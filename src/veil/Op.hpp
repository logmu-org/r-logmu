// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <array>
#include <cstdint>
#include <string_view>

namespace veil
{

// The veil operation set.
//
// Monikers are shared across base types where the operation means the same thing: `Add` covers both
// `double` and integer (datey / durationy) addition, and the result type follows from the argument
// types via `ResultRule`. Ops that are declared but not yet implemented are marked in the table;
// they are kept here so the shape is settled and filling one in is a local change.
enum class Op : uint16_t
{
    // Arithmetic -- unary
    Pos = 0,
    Neg,
    Abs,
    Recip,
    Sqrt,
    Exp,
    Log,
    Log10,
    Expm1,
    Log1p,
    Sin,
    Cos,
    Floor,
    Ceiling,
    Round,
    Trunc,
    Sign,
    Acc,

    // Arithmetic -- binary
    Add,
    Sub,
    Mul,
    RDiv,
    Pow,
    Min,
    Max,
    QScale,
    IDiv,
    IMod,
    Hypot,
    CopySign,

    // Arithmetic -- ternary
    Clamp,
    Fma,

    // Logical
    Not,
    And,
    Or,
    Xor,

    // Comparison
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,

    // NA test. What counts as NA is per base type: a double is NaN (veil draws no distinction
    // between R's NA_real_ and any other NaN), and a datey is any value below the representable
    // range. Bool and text cannot be NA in exp_data / val_data at all.
    IsNa,

    // Per-individual selection between two branches on a time-invariant condition. This is what
    // `if (cond) a else b` and `ifelse(test, yes, no)` become. Lowering turns it into the IL's
    // structured if/else; that the condition does not vary over time is what makes the choice once
    // per individual rather than once per slot, and is checked during tree rewriting.
    Select,

    // Conversion
    ToDouble,
    Broadcast,

    // Time vector construction
    VectorT,
    VectorDurn,
    VectorLogMu,

    // Finalisation -- collapse a time vector to a scalar
    Integrate,
    DiedValue,

    Count_
};

enum class OpCategory : uint8_t
{
    Arithmetic,
    Logical,
    Comparison,
    Selection,
    Conversion,
    TimeVector,
    Finalise,
};

// How the result type is derived from the argument types.
enum class ResultRule : uint8_t
{
    SameAsArgs,     // Arithmetic: the common argument type (Double or Int)
    SameAsBranches, // Selection: the common type of args 1 and 2; arg 0 is the condition
    AlwaysDouble,
    AlwaysBool,
    AlwaysVector,   // A `double` time vector
};

// WHICH TYPES AN OP ACCEPTS IS NOT RECORDED HERE. It was, as an `OpDomain` bitmask, and nothing ever
// read it: the real rules are relations between operands rather than membership per operand --
// `comparable` and `branchCommonType` in passAnnotateTypes, and the per-category storage checks in
// passLowerToBlock. A bitmask can state neither, so it sat inert while reading as a rule that was
// being enforced, which is how text ordering came to look approved. A type restriction belongs in
// the pass that can express it.
struct OpInfo final
{
    Op op;
    std::string_view moniker;
    uint8_t arity;
    OpCategory category;
    ResultRule resultRule;
    // Vectorisable ops accept a time vector in any argument; if any argument is a time vector then
    // so is the result. Non-vectorisable ops are scalar-only.
    bool vectorisable;
    // Read by the R binding, which refuses a declared-but-unbuilt op by name rather than letting it
    // fail later inside a pass.
    bool implemented;
};

namespace detail
{

constexpr std::array<OpInfo, static_cast<size_t>(Op::Count_)> opTable = {{
    // op             moniker        arity category                resultRule                vectorisable implemented
    {Op::Pos,         "pos",         1, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Neg,         "neg",         1, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Abs,         "abs",         1, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Recip,       "recip",       1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  false},
    {Op::Sqrt,        "sqrt",        1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Exp,         "exp",         1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Log,         "log",         1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Log10,       "log10",       1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Expm1,       "expm1",       1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Log1p,       "log1p",       1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Sin,         "sin",         1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Cos,         "cos",         1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Floor,       "floor",       1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Ceiling,     "ceiling",     1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Round,       "round",       1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Trunc,       "trunc",       1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Sign,        "sign",        1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Acc,         "acc",         1, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  false},

    {Op::Add,         "add",         2, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Sub,         "sub",         2, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Mul,         "mul",         2, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::RDiv,        "rdiv",        2, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Pow,         "pow",         2, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},
    {Op::Min,         "min",         2, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Max,         "max",         2, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::QScale,      "qscale",      2, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  false},
    {Op::IDiv,        "idiv",        2, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::IMod,        "imod",        2, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Hypot,       "hypot",       2, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  false},
    {Op::CopySign,    "copysign",    2, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  false},

    {Op::Clamp,       "clamp",       3, OpCategory::Arithmetic, ResultRule::SameAsArgs,   true,  true},
    {Op::Fma,         "fma",         3, OpCategory::Arithmetic, ResultRule::AlwaysDouble, true,  true},

    {Op::Not,         "not",         1, OpCategory::Logical,    ResultRule::AlwaysBool,   false, true},
    {Op::And,         "and",         2, OpCategory::Logical,    ResultRule::AlwaysBool,   false, true},
    {Op::Or,          "or",          2, OpCategory::Logical,    ResultRule::AlwaysBool,   false, true},
    {Op::Xor,         "xor",         2, OpCategory::Logical,    ResultRule::AlwaysBool,   false, true},

    {Op::Eq,          "EQ",          2, OpCategory::Comparison, ResultRule::AlwaysBool,   true,  true},
    {Op::Ne,          "NE",          2, OpCategory::Comparison, ResultRule::AlwaysBool,   true,  true},
    {Op::Lt,          "LT",          2, OpCategory::Comparison, ResultRule::AlwaysBool,   true,  true},
    {Op::Le,          "LE",          2, OpCategory::Comparison, ResultRule::AlwaysBool,   true,  true},
    {Op::Gt,          "GT",          2, OpCategory::Comparison, ResultRule::AlwaysBool,   true,  true},
    {Op::Ge,          "GE",          2, OpCategory::Comparison, ResultRule::AlwaysBool,   true,  true},

    {Op::IsNa,        "is_na",       1, OpCategory::Comparison, ResultRule::AlwaysBool,   true,  true},

    {Op::Select,      "select",      3, OpCategory::Selection,  ResultRule::SameAsBranches, true, true},

    {Op::ToDouble,    "double",      1, OpCategory::Conversion, ResultRule::AlwaysDouble, true,  true},
    {Op::Broadcast,   "vector",      1, OpCategory::Conversion, ResultRule::AlwaysVector, false, true},

    {Op::VectorT,     "vector_t",    0, OpCategory::TimeVector, ResultRule::AlwaysVector, false, true},
    {Op::VectorDurn,  "vector_durn", 1, OpCategory::TimeVector, ResultRule::AlwaysVector, false, true},
    {Op::VectorLogMu, "vector_log_mu", 2, OpCategory::TimeVector, ResultRule::AlwaysVector, false, true},

    {Op::Integrate,   "integrate",   1, OpCategory::Finalise,   ResultRule::AlwaysDouble, false, true},
    {Op::DiedValue,   "died_value",  1, OpCategory::Finalise,   ResultRule::AlwaysDouble, false, true},
}};

// The table is indexed by the enumerator, so every row must sit at its own index. Adding an op
// without adding its row -- or adding rows out of order -- fails the build rather than silently
// mislabelling an operation.
constexpr bool opTableIsWellFormed() noexcept
{
    for (size_t i = 0; i < opTable.size(); ++i)
    {
        if (static_cast<size_t>(opTable[i].op) != i) { return false; }
        if (opTable[i].arity > 3) { return false; }
        if (opTable[i].moniker.empty()) { return false; }
    }
    return true;
}

static_assert(opTable.size() == static_cast<size_t>(Op::Count_), "Every Op needs exactly one table row.");
static_assert(opTableIsWellFormed(), "The op table is misordered, or a row is malformed.");

} // namespace detail

constexpr const OpInfo& opInfo(Op op) noexcept
{
    return detail::opTable[static_cast<size_t>(op)];
}

constexpr std::string_view moniker(Op op) noexcept
{
    return opInfo(op).moniker;
}

// The six ordering/equality comparisons -- not `IsNa`, which shares their `OpCategory` but is unary
// and answers a different question.
constexpr bool isComparisonOp(Op op) noexcept
{
    return op == Op::Eq || op == Op::Ne || op == Op::Lt || op == Op::Le || op == Op::Gt || op == Op::Ge;
}

// The four ordering comparisons, as against `Eq` and `Ne`. A code records WHICH value, not where it
// sits in any sequence, so text is refused these -- in passAnnotateTypes, with passLowerToBlock
// restating it as a guard.
constexpr bool isOrderingComparison(Op op) noexcept
{
    return op == Op::Lt || op == Op::Le || op == Op::Gt || op == Op::Ge;
}

// `a OP b` read the other way round, so a tree-rewriting pass need only handle `field OP literal` and
// normalise `literal OP field` to it by reversing the operator (`Eq`/`Ne` are their own reverse).
// Defined only over `isComparisonOp`; callers must filter first.
constexpr Op reverseComparison(Op op) noexcept
{
    switch (op)
    {
        case Op::Lt: return Op::Gt;
        case Op::Le: return Op::Ge;
        case Op::Gt: return Op::Lt;
        case Op::Ge: return Op::Le;
        default: return op; // Eq, Ne.
    }
}

} // namespace veil
