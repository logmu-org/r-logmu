# logmu

R package providing high performance actuarial mortality experience analysis and modelling.

## Project outline

- This is a CRAN-targeted R package developed in RStudio
- C++ code lives in `src/` (but this is to be reviewed if we need to expose
functionality to other packages)
- This package uses the date and duration framework set out in the datey package
that I'm also developing and which has been submitted to CRAN

## Functionality for users

- In effect this is an implementation of the concepts set out in the vignettes 
now imported in the vignettes folder.

- Experience and valuation seriatim data -- these are identical other than
    - Experience data (`exp_data`) relates to an experience period and each individual has an
    E2R, comprising exposure start and end times and an indicator as to whether
    they died at the end of that exposure
    - Valuation data (`val_data`) is individual as at a single point in time

    A `datey` column called `birth` must always be present.

- The ability to define functions of the data using a pronoun convention where
    - `.i$FFF` accesses field `FFF`
    - `.t` accesses time (a `datey`)
    - `.b` is shorthand for `.i$birth`
    - `.x` is shorthand for `.t - .b`

- Define various S3 types that represent functions of data using a pronoun convention where `.i$FFF` accesses field `FFF` and `.t` accesses time (a `datey`)
    - `variable` -- a scalar function of `.i` and `.t`.
    - `mortality` -- a scalar function of `.i` and `.t` defined as $\log\mu(i,t)$.
    - `include` -- a `datey_interval` or `logical` function of `.i` only 
    analogous to base R's `subset`.

- The `variable` and `mortality` types should 'know':

    - Whether they depend on `.t`, `.b` or any `.i` fields other that `birth`.
    - If they depend on `.t`, their scale, i.e. the largest `durationy` numerical integration scalar that is safe to use.
    - A broad indication of the values they may take. For instance if we know a `variable` can only ever take the values 0 or 1 then
    we need calculate only Ew as opposed to Ew and E(w^2).

- For the initial release, the only `mortality` objects (i.e. not defined by pronoun expressions) I need are

    - Constant (and I'm not sure about this -- a pronoun expression can do this trivially)
    - Age-period tables defined by a 2D matrix of $\log\mu$ by annual age and year (with fixed `durationy`/`datey` offsets)

- For the subsequent releases, I may also need

    - `base_table` - $\log\mu$ by annual age (with a fixed `durationy` offsets) and a `datey` as at date
    - `projection` - $q$ or $\mu$-based age-period mortality projections defined by an age-period table of log projection factors
    - `variation` - can be layered in to a `mortality` to provide sub-annual mortality variation

- Provide a library system to be hooked into by dependent packages that use an R-like syntax that looks like this:

    ```R
    cmi::S4PMA # A CMI base table
    cmi::CMI_2025(0.015) # A CMI mortality projection that takes parameters
    ```
- Functions for users, all of which use NSE:

    - `aev(exp_data, mortality, include = NULL, weight = NULL)` calculates 
    $A=\text{A}w$, $E=\text{E}w$ and $V=\Omega\text{E}w^2$, 
    saving the result in the S3 type `aev`.
    $V$ is set to $\Omega E$ if the weight is 0 or 1 
    (i.e. an indicator, which is expressed in R as returning a logical).
    $\Omega$ is user-supplied overdispersion.

    - `aev_multi(exp_data, mortality, include = NULL, weight = NULL, sub_includes = NULL)`
    as same as `aev()` but returns a list of (named?) `aev`s. Provide 
    tidyverse-style A/E confidence interval charts with residuals.

    - Fit categorical values using statistics, which may be ordered or unordered.
    An example is optimal clustering of postcode-based socio-economic types.

    - `fit` fits proportional hazards mortality model (`model`) using maximum weighted log-likelihood.

    - Select between fitted mortality models using weighted AIC.

    - Scale standard tables as $q$ or $\mu$ scaling using PV equivalence to general mortalities (based on very simple and performant annuity calcs)

## Log-likelihood definitions

$$\begin{aligned}
\text{Var}\big(\text{A}w-\text{E}w\big) &= \Omega\,\mathbb{E}\,\text{E}w^2
\\
L^* &= \Omega^{-1}(\text{A}w\log\mu-\text{E}w)
\\
\mathbf{I} &= \Omega^{-1}\text{E}wXX^\text{T} = -L^*{}'' 
\\
\mathbf{J} &= \Omega^{-1}\text{E}w^2XX^\text{T} \mathrel{\hat=} \text{Var}(L^*{}') 
\\
L &= Z^{-1}L^*
\\
p &= Z^{-1}\text{tr}\big(\mathbf{J}\mathbf{I}^{-1}\big)
\\
L_\text{P} &= L(\hat\beta)-p
\\
\text{Var}\big(\hat\beta\big)&\mathrel{\hat=}\mathbf{I}^{-1}\mathbf{J}\mathbf{I}^{-1}
\end{aligned}$$

where

- $\Omega$ is user-supplied over-dispersion.

- $Z$ is $\text{E}w^2 / \text{E}w$ on a test mortality across the same experience data.
  Note that this is $V/(\Omega E)$ in the above AEV notation.

## Implementation

- High performance calculation of the functions for uses in the previous section
using SIMD vector code and multi-threading where appropriate

- Prior review of specialised calculations to sequence them into standard scalar and vectorised operations.
E.g. `aev(exp_data, include, mortality, weight)` would use NSE to examine the 
arguments to generate a set of operations that is then passed to a performant
calculation engine written in C++.

- The functions for users can be batched using the `batch()` function, which 
also uses NSE, and works something like this:

    Instead of
    
    ```R
    aev_m     <- aev(exp_data, mortality = mortality_m, include = include_m,     weight = amounts)
    aev_f_ret <- aev(exp_data, mortality = mortality_f, include = include_d_ret, weight = amounts)
    aev_f_dep <- aev(exp_data, mortality = mortality_f, include = include_d_dep, weight = amounts)
    ```

    the user writes
    
    ```R
    b <- batch(
      aev_m     = aev(exp_data, mortality = mortality_m, include = include_m,     weight = amounts),
      aev_f_ret = aev(exp_data, mortality = mortality_f, include = include_d_ret, weight = amounts),
      aev_f_dep = aev(exp_data, mortality = mortality_f, include = include_d_dep, weight = amounts)
    )
    ```

    and the results are run by the engine in one go and the results are in 
    `b$aev_m`, `b$aev_f_ret` etc.

## CRAN compliance

- Warn me if any changes will generate CRAN notes, warnings or errors on any
current CRAN platform when I submit the package via 
`devtools::check_win_devel()` or similar checks.
- Check Roxygen2 comments for CRAN fails such as non-ASCII characters and 
missing `@returns` statements.
- Avoid Undefined Behaviour (UB)
- Never use pointer casting or unions for type punning; use std::bit_cast
- Code must pass aggressive -Wstrict-aliasing=2 optimizations

## Standards

- R minimum version: **4.6.0**
- Required C++ standard: **C++20**
- All R and C++ files other than those sourced externally must have a standard 
disclaimer at the top. This excludes auto-generated files (e.g. `cpp11.R`) 
and vendored files, e.g. the `src/eve` or `src/spy` sub folders.

    The disclaimer for R files is:

    ```R
    # LogMu mortality experience analysis and model fitting
    #
    # This file is licensed to you under the Apache Licence 2.0.
    #
    # Copyright (c) Tim Gordon
    ```

    The disclaimer for C++ (including `.cpp` and `.hpp`) files is:

    ```cpp
    // LogMu mortality experience analysis and model fitting
    //
    // This file is licensed to you under the Apache Licence 2.0.
    //
    // Copyright (c) Tim Gordon
    ```

## Dependencies

- For R/C++ code interop use the `cpp11` package, i.e. `cpp11::` `doubles`, `logicals`, `sexp`, etc. and `[[cpp11::register]]`.
- Do *not* use `Rcpp`, `Rcpp::NumericVector`, or `Rcpp` macros.
- For multi-threading, use C++20's `std::jthread`.
- For high performance SIMD functions on vectors of double use vendored EVE SIMD library by Falcou et al (and its dependency SPY)

## Code style

- For R code follow Tidyverse naming conventions
- For C++ code in the `veil` sub-folder use
    - `PascalCase` for naming types and constants, and
    - `camelCase` for naming functions and function parameters
- For C++ code in general (other than vendored code, which *must* be left intact):
    - Use modern C++20 idioms
    - Use Allman bracket style
    - Check for const T&
    - Use range-based for loops, const auto&, and explicit single-argument constructors
    - Prefer wordy identifiers rather than abbreviations unless they are universal or keywords. For instance, prefer `literal` to `lit`.
    - Do *not* use indentation styles that rely on fixed width fonts. 
      Specifically do *not* align to opening delimiters (`AlignAfterOpenBracket: Align` in CLang, `Align contents to opening parenthesis` Visual Studio etc).
      Instead used hanging / block indent with a fixed indentation level regardless of the function name's length.
      For example, do *not* do this:

      ```cpp
      int some_function(int a,
                        int b,
                        int c)
      ```
      
      but *do* do this:

      ```cpp
      int some_function(
        int a,
        int b,
        int c)

## Language style

Act as a direct, clear human writer.

Always use British English and punctuation. In particular:
- When using a dash as a strong comma, colon, or parenthesis to create an emphatic pause, break in thought, or aside within a sentence, always use en dash (–) with spaces either size. Never use em dash (—).
- Use single quotation marks before double quotation marks.
- Never use 'off of'.
- Always use 'named after' instead of 'named for' (unless you are giving the reason for the naming).

Do not use AI-tell phrasing:

- Do not use telegraphic appositive headlines or fragments ("One construct, three uses", "Same core, two bindings").
- Do not use rule-of-three parallelism for rhythm rather than content.
- Do not use "it isn't X, it's Y" reframes or em or en-dash reframes deployed for punch.
- Do introduce lists and points with plain lead-in clauses.
- Do not use any of these words: delve, realm, harness, unlock, tapestry, beacon, testament, symphony, navigate, journey, furthermore, moreover, ultimately, pivotal, broader, landscape.
- Do not use overly complex punctuation.
- Do not open with introductory filler or close with a summary.
- Vary your sentence length drastically -- make some punchy and one sentence long.
- Do not explain a concept and then immediately summarise it.

All that said, README.Rmd is in effect a marketing page and is permitted to be salesy.

## Things Claude should not do

- Do not amend this file `CLAUDE.md`.
- Do not request permission to amend this file `CLAUDE.md`.
- Do not amend any code without my permission.
- Do not amend the vendored code in `src/eve` or `src/spy` without my express and clear permission.
- Do not download any software.
- Do not commit any changes without my permission.
- Do not offer to commit changes yourself (although you can suggest that I do commit changes at sensible stopping points).
- Do not add Claude credits e.g. "Co-Authored-By: Claude" in code comments or in git commit summaries or descriptions.
