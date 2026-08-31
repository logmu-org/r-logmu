// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <stop_token>
#include <thread>
#include <vector>

namespace veil
{

// Running a fixed set of independent tasks across threads.
//
// THIS CANNOT MOVE AN ANSWER, AND THAT IS THE WHOLE POINT. The caller has already divided its work
// into numbered tasks and has somewhere to put each task's result, so all that happens here is that
// a different thread may walk a given task. Nothing about the schedule reaches a total, because the
// totals are folded by task index afterwards -- see RecordChunk.hpp for why the division is arranged
// that way. The test this earns is the strongest one available: the same data at any thread count
// must give bit-identical results.
//
// THE SHAPE IS THE EASIEST A POOL CAN HAVE, so it stays small. The task set is known up front, tasks
// are independent, none creates more work, none blocks on I/O, and there is no nesting. That rules
// out everything a general pool needs -- work stealing, futures, a task graph, a queue per worker --
// and leaves one atomic counter that workers increment.
//
// C++ HAS NO STANDARD THREAD POOL, and the near misses are traps. `std::async` may or may not pool
// depending on the implementation, and libstdc++ -- which is what Rtools gives us -- spawns a fresh
// thread every call. The C++17 parallel algorithms need TBB under libstdc++. `std::execution` is
// C++26. Thirty lines of our own is smaller than any of those dependencies.

// How many worker threads to use for a given request.
//
// A COUNT RATHER THAN A FLAG. A flag cannot say "two" for a CRAN check or "four" on a shared server,
// and both of those are real requests. Exposing a count at all is only safe BECAUSE the answer does
// not depend on it: turning threading off to check something must not move a digit, and it does not.
//
// 1 means run on the calling thread and spawn nothing. 0 means as many as the machine reports, which
// is an explicit opt-in and never a default. Anything else is taken literally.
//
// `hardware_concurrency` IS A HINT AND MAY RETURN 0, so it is always guarded. It also counts logical
// processors, which for floating-point-bound work over-states the useful parallelism, and it
// typically ignores affinity masks and container CPU limits. There is no standard way to ask for
// physical cores, so this is what there is.
inline size_t resolveThreadCount(size_t requested) noexcept
{
    if (requested != 0) { return requested; }

    const unsigned int reported = std::thread::hardware_concurrency();
    return reported == 0 ? 1 : static_cast<size_t>(reported);
}

// The most threads a caller may take without asking.
//
// CRAN CHECKS ON TWO CORES, and a package that helps itself to more during a check is a note at
// best. The cap belongs here rather than in the R binding so that every front end inherits it.
constexpr size_t DefaultThreadCount = 2;

// What a parallel run produced. The results themselves live wherever the task body put them; this is
// only how the run ended.
struct ParallelRunOutcome final
{
    // False when the host asked to stop before every task had run, in which case the task results are
    // incomplete and the caller should raise rather than report them.
    bool completed = true;
};

// Runs `taskCount` tasks, calling `body(taskIndex)` once for each, and returns once every task has
// run or the host has asked to stop.
//
// THE CALLING THREAD PULLS TASKS TOO. Spawning n workers and then blocking would idle a core, and on
// a two-core CRAN check that is half the throughput. It also means the calling thread is the one
// thread whose position in the work is known, which is what makes the host poll below safe.
//
// `hostPoll` IS CALLED ONLY BY THE CALLING THREAD, and only between its own tasks. It answers true
// when the host wants the run stopped. This is how a host interrupt reaches a core that must not know
// what a host is: R's `R_CheckUserInterrupt` is main-thread-only and longjmps, so the R binding polls
// for a pending interrupt WITHOUT jumping, answers true, and raises the R error only after this
// function has returned and every thread has been joined. Nothing R-specific enters the core.
// It may be empty, in which case nothing is polled and the run cannot be stopped early.
//
// EXCEPTIONS ARE CAUGHT PER TASK, not per worker, and the FIRST BY TASK INDEX is rethrown. One
// escaping a thread function calls `std::terminate`, so catching is not optional -- but catching per
// worker would make which error surfaces depend on which worker happened to take which task, and an
// error message that varies with the schedule is the same disease as a total that does. A task index
// is the caller's own numbering, so the error is as reproducible as the answer would have been.
//
// Interrupt latency is one task: a worker checks the stop token between tasks, never inside one.
template <typename TaskBody>
inline ParallelRunOutcome runInParallel(
    size_t taskCount,
    size_t threadCount,
    const std::function<bool()>& hostPoll,
    const TaskBody& body)
{
    ParallelRunOutcome outcome;
    if (taskCount == 0) { return outcome; }

    const size_t workers = resolveThreadCount(threadCount);

    // NOT WORTH SPAWNING FOR. A single task cannot be shared out, and a single worker is the calling
    // thread by definition. Both run inline, with no thread created and no atomic touched, which also
    // means a threading bug cannot reach a caller who turned threading off.
    if (workers <= 1 || taskCount == 1)
    {
        for (size_t task = 0; task < taskCount; ++task)
        {
            if (hostPoll && hostPoll())
            {
                outcome.completed = false;
                return outcome;
            }
            body(task);
        }
        return outcome;
    }

    std::atomic<size_t> nextTask{0};
    std::stop_source stopSource;

    // One slot per task rather than per worker, so that the error that surfaces is the lowest-
    // numbered failing task whichever worker happened to run it.
    std::vector<std::exception_ptr> failures(taskCount);

    const std::stop_token stopToken = stopSource.get_token();
    const auto runTasks = [&]()
    {
        while (!stopToken.stop_requested())
        {
            const size_t task = nextTask.fetch_add(1, std::memory_order_relaxed);
            if (task >= taskCount) { return; }

            try
            {
                body(task);
            }
            catch (...)
            {
                failures[task] = std::current_exception();
                stopSource.request_stop();
                return;
            }
        }
    };

    {
        // SPAWNED PER CALL, DELIBERATELY. Eight threads cost a couple of hundred microseconds to
        // create, which is nothing beside a real calculation. A pool kept alive between calls would
        // buy sleeping threads in every R session, a lifetime to manage, and a genuine deadlock
        // hazard: forking a process that has threads is a known trap, and `parallel::mclapply` forks.
        //
        // `std::jthread` joins in its destructor, so this scope is the join. The workers are created
        // before the calling thread starts pulling so that they are not waiting on it.
        std::vector<std::jthread> pool;
        pool.reserve(workers - 1);
        for (size_t worker = 0; worker + 1 < workers; ++worker) { pool.emplace_back(runTasks); }

        // The calling thread's own share, interleaved with the host poll.
        while (!stopToken.stop_requested())
        {
            if (hostPoll && hostPoll())
            {
                outcome.completed = false;
                stopSource.request_stop();
                break;
            }

            const size_t task = nextTask.fetch_add(1, std::memory_order_relaxed);
            if (task >= taskCount) { break; }

            try
            {
                body(task);
            }
            catch (...)
            {
                failures[task] = std::current_exception();
                stopSource.request_stop();
                break;
            }
        }
    }

    // Every worker has been joined by here, so the failure slots are stable and nothing is still
    // running that could touch the caller's results.
    for (const std::exception_ptr& failure : failures)
    {
        if (failure) { std::rethrow_exception(failure); }
    }

    return outcome;
}

} // namespace veil
