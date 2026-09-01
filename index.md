# logmu

If you are setting or reviewing base mortality for defined benefit
pension, annuity or life assurance liabilities then **logmu** is
designed for you[^1], with state of the art tooling, including support
for postcode-based models, and performance at scale.

The **logmu** feature set:

- **Flexible mortality framework** incorporating time-based covariates
  and arbitrary proportional hazards models.

- **A/E analytics** including confidence intervals and residuals, with
  visualisation.

- **Fit proportional hazards models** by maximum likelihood, select
  between them by AIC, and cluster categorical covariates — all in a
  single pipeline.

- **Customisation** — weighting (e.g. amounts vs lives), statistical
  relevance (e.g. treat older data as less reliable) and time-based
  inclusion criteria (e.g. sub-setting experience by age) are built into
  every function.

- **Blazingly fast** – **logmu** takes full advantage of SIMD
  vectorisation and multi-threading. A/E for **XX** thousand records
  takes under Y seconds.

The principles underlying **logmu** have been applied by the author to
longevity transactions totalling c.£100Bn across the UK, US and
Netherlands.

A typical A/E analysis — data to charts — takes under 20 lines of R:

``` r

# exp_data: experience dataset with fields sex (character) and pension (double)
S4PMA <- mortality_table(...)
S4PFA_110 <- mortality_table(...)

mortality <- mortality(if (.i$sex == "male") mort_m else mort_f)
includes  <- aged(60, by = 5, to = 100) | periods(2010, by = 5, to = 2025)

by_lives   <- aev(exp_data, mortality, includes, weight = .i$pension > 0)
by_amounts <- aev(exp_data, mortality, includes, weight = .i$pension)

plot(by_lives)
plot(by_amounts)
```

[Get started](https://r-logmu.logmu.org/articles/A000-get-started.html)

[^1]: For avoidance of doubt, in case you’ve been led here by searching
    on common mortality keywords, **logmu** is not (currently) for
    mortality *trend* analysis – the techniques and data make it a
    distinct discipline, or *medical mortality statistics* – try
    [survival](https://CRAN.R-project.org/package=survival).
