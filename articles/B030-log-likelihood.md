# Log-likelihood

This article was originally published on
[timgord.com](https://timgord.com/2025-08/mortality-log-likelihood/) on
2025-08-12.

I think it’s a shame that the ‘log’ in ‘log-likelihood’ is so often
presented as a technical convenience or a device for avoiding numerical
under/overflow. Yes, it is definitely both of these things, but it is
much more fundamental.

**Expected *log*-probability, i.e. *entropy*, lies at the heart of
information theory.** And the concept of entropy itself is pervasive,
having extended beyond thermodynamics, its original home, into quantum
physics and general relativity, as well as information theory.

So, without further ado, let’s define *log*-likelihood for mortality
experience data.

## Definition

Although we can approach log-likelihood from an all-of-time point of
view, I think it’s instructive to start at the opposite end of the
scale, with infinitesimals[^1].

Over an infinitesimal time period \text{d}t, the probability of survival
is

p\approx\exp\Big\\\\-\\\mu(i,t)\text{d}t\Big\\

There are two possible outcomes, which I’ll represent using \delta,
which is 1 if the individual died and 0 otherwise. The likelihood over
this infinitesimal time period is then

p^{(1-\delta)}(1-p)^\delta\approx\exp\Big\\\\-\\(1-\delta)\mu(i,t)\text{d}t\Big\\
\Big\\\mu(i,t)\text{d}t\Big\\^\delta

and hence the *log*-likelihood is

\delta\log\Big\\\mu(i,t)\text{d}t\Big\\ - (1-\delta)\mu(i,t)\text{d}t

How can we take the log of the infinitesimal \text{d}t? The answer is
that, even though they are often not presented this way, log-likelihoods
are always *relative*, i.e. it is only the *difference* between two
log-likelihoods that matters. And, when we subtract another
log-likelihood *for the same experience data*[^2], the \text{d}ts inside
the log term will cancel.

Given our assumption of independence by individual and time period, we
can sum these individual log-likelihoods. When we do this, \delta is, by
definition zero everywhere except the very end and so the (1-\delta) can
be treated as 1, which leaves us with[^3]

\delta\log \mu(i,t) - \mu(i,t)\text{d}t\tag{3}

But we already know how to add up these terms using the [\text{A} and
\text{E}
operators](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Def-AE-ops)
defined in the [first
article](https://r-logmu.logmu.org/articles/B010-measures-matter.md) in
this series. So the log-likelihood is

L=\text{A}w\log\mu-\text{E}w\tag{4}

where w is a
[variable](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Def-variable)
weighting of the log-likelihood infinitesimals.

\text{A} and \text{E} are
[measures](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Def-measure),
so L is a measure[^4] too, which means that, just like \text{A} and
\text{E}, we can partition the data any way we like – we’ll still end up
the same result
([Insight 2](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Insight2)).

## Weighting

Log-likelihood is literally the log of a probability so the inclusion of
a weight w implies that we’re using probabilities to the power of w,
which is worthy of comment.

1.  If w\in\\0,1\\ then this is simply equivalent to excluding or
    including data. I’ll call this ‘lives-weighted’.

2.  If 0\le w \le 1 then w can be interpreted as *relevance* (also known
    as ‘reliability’ or ‘importance’).

3.  The general case, w\ge0, is sometimes described as ‘*ad hoc*’ or
    ‘pragmatic’, or even illegitimate(!).

I’ll have more to say about weighting the log-likelihood in due course.
For the time being, let’s leave our options open by assuming w\ge0,
i.e. 3 above.

#### Insight 7. Log-likelihood can be defined directly in terms of the \text{A} and \text{E} operators

The log-likelihood written in terms of the \text{A} and \text{E}
operators is

L=\text{A}w\log\mu-\text{E}w

where w\ge0 is the weight variable.

(This is *before* allowing for overdispersion.)

#### Next article: [*Proportional hazards*](https://r-logmu.logmu.org/articles/B040-proportional-hazards.md)

Equation (4) is doing a lot of heavy lifting with admirable concision.
In the [next
article](https://r-logmu.logmu.org/articles/B040-proportional-hazards.md),
I’ll show how it leads directly to one of the most useful tools in the
mortality modelling armoury.

[^1]: My mental model of log-likelihood is as a sum of infinitesimals.
    It is all heuristics though – I am not making a pretence of
    mathematical rigour.

[^2]: A cardinal rule of log-likelihoods is that they are *comparable
    only if calculated on exactly the same data*.

[^3]: Taking the log of the [dimensional
    quantity](https://en.wikipedia.org/wiki/Dimensionless_quantity) \mu
    may also irk you – it does me (a little) – but units also cancel
    when two log-likelihoods are compared.

[^4]: Technically it’s a [signed
    measure](https://en.wikipedia.org/wiki/Signed_measure).
