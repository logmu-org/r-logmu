# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# THE CHART SIGNS ITS OWN THEME so that it can tell when the theme has been
# replaced. `theme_aev()` sets `logmu.signature`, and `facet_aev()` reads it back
# off the finished theme when the chart is drawn.
#
# WHY A THEME ELEMENT AND NOT AN ATTRIBUTE. The signature has to survive the one
# thing that is not a problem -- `+ theme(...)`, which MERGES and is how anybody
# customises a chart -- while being lost to the thing that is: a COMPLETE theme
# such as `ggplot2::theme_minimal()`, which replaces wholesale. A registered
# element does exactly that, measured on ggplot2 4.0.3:
#
#     theme_aev()                            TRUE
#     + theme(text = element_text(size=18))  TRUE   <- merged, still ours
#     + theme_minimal()                      NULL   <- replaced
#     + theme_bw()                           NULL
#
# An attribute would not: `add_theme()` builds a new list and does not carry it.
#
# `register_theme_elements()` is ggplot2's documented extension point and writes
# into its global element tree, which is why this belongs in `.onLoad` and why
# `ggplot2::reset_theme_settings()` undoes it.
.onLoad <- function(libname, pkgname) {

  ggplot2::register_theme_elements(
    logmu.signature = FALSE,
    element_tree = list(logmu.signature = ggplot2::el_def("logical"))
  )

  invisible(NULL)
}
