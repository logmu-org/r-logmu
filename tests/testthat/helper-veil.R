# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The permitted integration time scales, in clicks, named the way these tests
# talk about them.
#
# ROUTED THROUGH `time_scale_clicks()` rather than written out as numbers, so a
# test validates against the engine's own list. A literal 133590 here would be a
# second place for the permitted set to live, and the failure it hides is the
# one the C++ guard exists to catch: the entry points take these positionally,
# so a wrong value is a plausible answer rather than an error.
# QUALIFIED, unlike everything in the test files themselves. `test_check()` runs
# those inside the package namespace, so internals are visible there without a
# prefix -- but a helper is sourced BEFORE that happens and does not see them.
year_scale <- logmu:::time_scale_clicks(1)
quarter_scale <- logmu:::time_scale_clicks(1 / 4)
month_scale <- logmu:::time_scale_clicks(1 / 12)
# The finest permitted scale, a fifth of a month. No common name, hence this one.
fine_scale <- logmu:::time_scale_clicks(1 / 60)

# The overdispersion the engine tests run at.
#
# ONE, DELIBERATELY AND EVERYWHERE IN THESE FILES. V then comes back as the
# plain second moment Ew^2, which is what every analytic oracle here is written
# against, and it is what keeps the `V = E` assertions for an indicator weight
# correct. A test that wants to see overdispersion do something must say so.
no_overdispersion <- 1
