# Parse a pronoun expression into an AST

Parses an R expression written with the data pronouns `.i$field`, `.t`,
`.b` and `.x` into a small abstract syntax tree. Either a bare
expression or a `~` formula may be supplied; constant sub-expressions
are evaluated in the calling environment (or the formula's environment)
and folded.

## Usage

``` r
pronoun_expressions(expression)
```

## Arguments

- expression:

  A pronoun expression, optionally written as a `~` formula.

## Value

An `it_node`, the root of the parsed expression tree.

## Examples

``` r
pronoun_expressions(.i$pension > 0)
#> call `>`
#>   field .i$pension
#>   lit 0
pronoun_expressions(~ .t - .b)
#> call `-`
#>   time .t
#>   field .i$birth
```
