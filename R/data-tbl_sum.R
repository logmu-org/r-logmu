# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#### exp_data ####

#' @importFrom pillar tbl_sum
#' @export
tbl_sum.exp_data <- function(x) {

  start <- exp_start(x)
  end <- exp_end(x)
  duration <- end - start

  c(
    #"Exp data" = pillar::dim_desc(x),
    "Exp data" = paste0(pillar::dim_desc(x), " (", duration, ")"),
    #"\u00A0\u00A0Period" = paste0("[", format(exp_start(x)), ", ", format(exp_end(x)), ")")
    "\u00A0\u00A0Start" = format(start),
    "\u00A0\u00A0End" = format(end)
  )
}


#### val_data ####

#' @importFrom pillar tbl_sum
#' @export
tbl_sum.val_data <- function(x) {
  c(
    "Val data" = pillar::dim_desc(x),
    "\u00A0\u00A0As at" = format(as_at(x))
  )
}
