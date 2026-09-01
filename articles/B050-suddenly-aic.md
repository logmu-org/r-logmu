# Suddenly AIC

This article was originally published on
[timgord.com](https://timgord.com/2025-08/mortality-suddenly-aic/) on
2025-08-26.

In this article, I’m going to look at choosing between mortality models
using Hirotugu Akaike’s information criterion, the AIC.

I’m going to run through – at a very high level – the rationale behind
the AIC and its construction because

- the standard result looks so trivial that people sometimes assumes
  it’s an arbitrary convention, and
- I’m going to generalise it (a little).

I [previously
wrote](https://r-logmu.logmu.org/articles/B030-log-likelihood.html#Ref-log)
that it’s a shame that the ‘log’ in ‘log-likelihood’ is not always
recognised as being much more fundamental than a mere technical
convenience or device for avoiding numerical under/overflow. The AIC is
a case in point – indeed Akaike himself wrote that ‘the log-likelihood
is essentially a more natural quantity than the simple likelihood’[^1].

## Akaike’s big idea

Selecting parameters by maximising the log-likelihood is an effective
tool for calibrating a single model. But we can’t choose between
different models by comparing maximum likelihoods because models with
more parameters can achieve a better fit to the data and so the
comparison is not fair.

#### Further reading

Kenneth Burnham and David Anderson’s [book on model
selection](https://doi.org/10.1007/b97636)[^2] is in my view the best
reference on *practical application* of the AIC. Also see (both online)
[AIC vs
BIC](http://www.sortie-nd.org/lme/Statistical%20Papers/Burnham_and_Anderson_2004_Multimodel_Inference.pdf)
and [AIC myths and
misunderstandings](https://sites.warnercnr.colostate.edu/kenburnham/wp-content/uploads/sites/25/2016/08/AIC-Myths-and-Misunderstandings.pdf).

Rob Hyndman provides an excellent overview of [facts and fallacies of
the AIC](https://robjhyndman.com/hyndsight/aic/).

Akaike’s[^3] insight was

- [*relative
  entropy*](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence)
  (also known as the *Kullback-Leibler divergence*) can be used to
  compare different models,

- *maximum likelihood is a biased estimate* of the variation in relative
  entropy compared with reality[^4],

- provided a model is reasonably good, *we can adjust for that bias with
  a simple penalty*, and

- the resulting *penalised maximum log-likelihood*, L\_\text{P}, can be
  used to compare models *regardless of the number of parameters they
  have*.

#### Too many twos

The AIC is defined as âˆ’2 times L\_\text{P}, for consistency with
regression.

Unfortunately the âˆ’2 factor

- complicates the definition,
- obscures the AIC as an unbiased estimate of [relative
  entropy](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence),
  and
- leads to spurious additional factors of 2 or Â½ in applications[^5].

So I’ll use the penalised log-likelihood, L\_\text{P}, itself
i.e. *without the âˆ’2*.

## Penalised log-likelihood

Although I’ll continue to refer to ‘the AIC’, I’ll express it in terms
of penalised log-likelihood, which is equivalent but cleaner (see box
out).

I’m not going to justify or derive the AIC here – if you want further
information then see Burnham & Anderson[^6] – but we will need a few key
results.

The penalised log-likelihood that falls out of the definition of
[relative
entropy](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence)
is

L\_\text{P} \mathrel{\hat=} \mathbb{E}\_1\mathbb{E}\_2\\L\tag{10}

where

- \mathbb{E} indicates expectation according to reality, which occurs
  once in the definition of relative entropy itself and a second time
  because we want our estimate to be unbiased, and

- L is log-likelihood (weighted by the
  [variable](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Def-variable)
  w\ge0), which we [previously
  defined](https://r-logmu.logmu.org/articles/B030-log-likelihood.html#Def-log-likelihood)
  using the [\text{A} and \text{E}
  operators](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Def-AE-ops)
  as

  L=\text{A}w\log\mu-\text{E}w\tag{4}

With some assumptions, we can estimate this as

L\_\text{P} = L(\hat\beta)-\operatorname{tr}\Big\\\mathbf{I} \\
\text{Var}\big(\hat\beta\big)\Big\\\tag{11}

where \beta parametrises the [mortality
model](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Def-mortality)
under consideration, \hat\beta maximises the log-likelihood, and all of
the following are \dim(\beta)^2 square matrices:

\begin{aligned}\text{Var}\big(\hat\beta\big)&\approx\mathbf{I}^{-1}\mathbf{J}\mathbf{I}^{-1}\\\rule{0pt}{3.5ex}\mathbf{I}&=-L''(\hat\beta)\\\rule{0pt}{3.5ex}\mathbf{J}&\mathrel{\hat=}\mathbb{E}\\L'(\beta)L'(\beta)^\text{T}\end{aligned}

with ' indicating \partial /\partial \beta and ^\text{T} indicating
transpose.

## The pay off

The \mathbf{I} matrix and its inverse cancel in equation (11), resulting
in

L\_\text{P}=L(\hat\beta)-\text{tr}(\mathbf{J}\mathbf{I}^{-1})\tag{12}

where \text{tr} is the [trace
operator](https://en.wikipedia.org/wiki/Trace_(linear_algebra)).

The usual next step is to note that, in the lives-weighted case,
i.e. w\in\\0,1\\, \mathbf{I}\approx\mathbf{J}, which results in the
conventional form[^7] of the AIC:

L\_\text{P}=L(\hat\beta)-\dim(\beta)\tag{lives-weighted only}

But we can still obtain a useful result *in the ad hoc weighted case* if
we assume that our mortality model is [proportional
hazards](https://r-logmu.logmu.org/articles/B040-proportional-hazards.html#the-proportional-hazards-model),
i.e.

\mu(\beta) = \mu^\text{ref}\exp\Big(\beta^\text{T}X\Big)\tag{9}

From
[before](https://r-logmu.logmu.org/articles/B040-proportional-hazards.html#Def-L-dash-dash),
we have

\mathbf{I}=-L''(\hat\beta)=\text{E}wXX^\text{T}

and, using the expected and variance results from the [A over E
article](https://r-logmu.logmu.org/articles/B020-a-over-e.md), we have

\begin{aligned}\mathbb{E}\\L'(\beta)L'(\beta)^\text{T}&=\mathbb{E}\big(\text{A}wX-\text{E}wX\big)\big(\text{A}wX-\text{E}wX\big)^\text{T}\\\rule{0pt}{3ex}&\approx\text{E}w^2XX^\text{T}\end{aligned}

resulting in our estimator for \mathbb{E}\\L'(\beta)L'(\beta)^\text{T}
being

\mathbf{J}=\text{E}w^2XX^\text{T}\tag{13}

and hence

\text{Var}\big(\hat\beta\big)\mathrel{\hat=}
\big(\text{E}wXX^\text{T}\big)^{-1}\big(\text{E}w^2XX^\text{T}\big)\big(\text{E}wXX^\text{T}\big)^{-1}\tag{14}

and

L\_\text{P}=L(\hat\beta)-\text{tr}\left(\frac{\text{E}w^2XX^\text{T}}{\text{E}wXX^\text{T}}\right)\tag{15}

#### Insight 9. An estimate of the variance of the fitted parameters for a proportional hazards mortality model is available in closed form for any *ad hoc* log-likelihood weight

\text{Var}\big(\hat\beta\big)\mathrel{\hat=}
\mathbf{I}^{-1}\mathbf{J}\mathbf{I}^{-1}

where \hat\beta is the maximum likelihood estimator of the covariate
weights, X is the vector of covariates, w\ge0 is the log-likelihood
weight, \mathbf{I}=\text{E}wXX^\text{T} and
\mathbf{J}=\text{E}w^2XX^\text{T}.

(This is *before* allowing for overdispersion.)

Caveat: w is an *ad hoc* reallocation of log-likelihood; it is *not*
[relevance](https://timgord.com/2025-10/mortality-good-things-come-to-those-who-weight-i/#3-defining-and-incorporating-data-relevance).
For the version of this insight that *does* take account of relevance,
see
[Insight 17](https://timgord.com/2025-11/mortality-good-things-come-to-those-who-weight-iii/#Insight17).

#### Insight 10. A penalised log-likelihood for a proportional hazards mortality model is available in closed form for any *ad hoc* log-likelihood weight

L\_\text{P}= L(\hat\beta)-\text{tr}\big(\mathbf{J}\mathbf{I}^{-1}\big)

where \hat\beta is the maximum likelihood estimator of the covariate
weights, X is the vector of covariates, L is the log-likelihood, w\ge0
is the log-likelihood weight, \mathbf{I}=\text{E}wXX^\text{T} and
\mathbf{J}=\text{E}w^2XX^\text{T}.

Caveat: w is an *ad hoc* reallocation of log-likelihood; it is *not*
[relevance](https://timgord.com/2025-10/mortality-good-things-come-to-those-who-weight-i/#3-defining-and-incorporating-data-relevance).
For the version of this insight that *does* take account of relevance,
see
[Insight 17](https://timgord.com/2025-11/mortality-good-things-come-to-those-who-weight-iii/#Insight17).

## Just weight a moment

*Provided we use a [proportional
hazards](https://r-logmu.logmu.org/articles/B040-proportional-hazards.html#the-proportional-hazards-model)
model*, we have formulas for penalised log-likelihood and the variance
of \hat\beta in the pragmatic *ad hoc* weighted case.

In the lives-weighted case, w\in\\0,1\\, then w^2=w and hence equation
(15) collapses to the lives-weighted version. Otherwise, weighting the
experience data will increase concentration and so the estimated
variance of \hat\beta will typically be greater than the lives-weighted
version. If that reflects the impact of \hat\beta on the liabilities
then this is intuitively reasonable.

There are issues though:

- The validity of equation (13) depends on the interpretation of the
  weight, w. If it is lives-weighted, i.e. w\in\\0,1\\, then everything
  is trivially ok. And if w is an *ad hoc* reallocation of
  log-likelihood then this seems ok too, albeit in a pragmatic hand-wavy
  sense[^8]. But equation (13) cannot be correct if w is interpreted as,
  for instance,
  [relevance](https://timgord.com/2025-10/mortality-good-things-come-to-those-who-weight-i/#3-defining-and-incorporating-data-relevance)
  – if it were then, say, halving it should double
  \text{Var}\big(\hat\beta\big), but it doesn’t. Hence the caveats.

- Being able to rank models by L\_\text{P} gets us only halfway because
  we also need to understand *whether differences in L\_\text{P} are
  significant*. In the lives-weighted and (after we’ve corrected for the
  point above) the relevance-weighted cases, this is straightforward – a
  difference in the penalised log-likelihood of 1 is significant because
  it is equivalent to adding one more parameter. But, for an *ad hoc*
  weight, using the formula as presented here does not automatically
  provide a measure of significance and so we need another way in[^9].

On balance, I think equation (15) is a useful practical tool, but we can
do better and put matters on a sounder footing. So I’ll revisit
weighting and relevance in a few articles’ time.

#### Next article: [*Overdispersion and quasi-log-likelihood*](https://r-logmu.logmu.org/articles/B060-overdispersion.md)

I used variance in the above *without allowing for overdispersion*. This
is [*always wrong in mortality
analysis*](https://r-logmu.logmu.org/articles/B010-measures-matter.html#Insight1),
so let’s tackle it in the [next
article](https://r-logmu.logmu.org/articles/B060-overdispersion.md).

[^1]: Akaike, H. (1973), “Information theory and an extension of the
    maximum likelihood principle”, in Petrov, B. N.; CsÃ¡ki, F. (eds.),
    *2nd International Symposium on Information Theory, Tsahkadsor,
    Armenia, USSR, September 2-8, 1971*, Budapest: AkadÃ©miai KiadÃ³,
    pp. 267–281. Republished in Kotz, S.; Johnson, N. L., eds. (1992),
    *Breakthroughs in Statistics*, vol. I, Springer-Verlag, pp. 610–624.

[^2]: Burnham, K. P.; Anderson, D. R. (2002), *Model Selection and
    Multimodel Inference: A practical information-theoretic approach*
    (2nd ed.), Springer-Verlag.
    <doi:%5B10.1007/b97636>\](<https://doi.org/10.1007/b97636>),
    ISBN-13: 9780387953649.

[^3]: Akaike, H. (1973), “Information theory and an extension of the
    maximum likelihood principle”, in Petrov, B. N.; CsÃ¡ki, F. (eds.),
    *2nd International Symposium on Information Theory, Tsahkadsor,
    Armenia, USSR, September 2-8, 1971*, Budapest: AkadÃ©miai KiadÃ³,
    pp. 267–281. Republished in Kotz, S.; Johnson, N. L., eds. (1992),
    *Breakthroughs in Statistics*, vol. I, Springer-Verlag, pp. 610–624.

[^4]: This does not imply that there is a ‘true model’ – models are by
    definition simplifications of reality.

[^5]: Examples:

    - The natural definition and standard guidance on what constitutes a
      significant difference in AIC is one parameter’s worth. This is 1
      for penalised log-likelihood, but a non-intuitive 2 for the
      conventional AIC definition.

    - Akaike weights, \exp(-\tfrac{1}{2}\text{AIC}), are pure clumsiness
      compared with \exp(L\_\text{P}).

[^6]: Burnham, K. P.; Anderson, D. R. (2002), *Model Selection and
    Multimodel Inference: A practical information-theoretic approach*
    (2nd ed.), Springer-Verlag.
    <doi:%5B10.1007/b97636>\](<https://doi.org/10.1007/b97636>),
    ISBN-13: 9780387953649.

[^7]: Barring a factor of âˆ’2.

[^8]: There is *non*-hand-wavy version to come in a few articles’ time.

[^9]: One option is to add up all the penalties of all candidate models
    and divide that the by the total number of parameters being fitted
    to obtain an equivalent to ‘one more parameter’.
