# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#### Where the numbers in this file come from ####
#
# The chart is a reimplementation of a prototype, and every constant below is
# taken from one of two places rather than chosen here:
#
# * `Chart prototypes - slim.pptx` in the `logmu-docs` repository, slides 6 and
#   7, which are marked "USE ME!!". They fix the scale, the marker, the layout
#   rules and the house rules for displayed text.
# * `AEChart.cs` and `ChartPalette.cs` in the `logmu-dev` repository, which are
#   the renderer the slides were produced from, and which give the geometry to
#   the pixel and the palette to the byte.
#
# WHERE THE TWO DISAGREE THE SLIDES WIN, because they are the design and the C#
# is one rendering of it. The C# maths is out of date in places and is not used
# here at all -- everything plotted is read off an `aev` through the accessors.

#### Palette ####
#
# `ChartPalette.Default`, the light theme. The three marker greys are the ones
# both sources agree on exactly; the rest are the C#'s, and reconciling them
# with the slides is a later job. `axis_line` is the C#'s `AxisLineColor`, which
# is `GridLineColor` put through `Darkened()` -- squaring L* in CIELAB -- and so
# is stored here as the colour that comes out rather than as the recipe.
aev_palette <- list(
  marker_weak     = "#B2B2B2",
  marker_medium   = "#969696",
  marker_strong   = "#4D4D4D",
  # THE CAP IS DARKER THAN THE STRONG GREY, and has its own entry so that it can
  # be. It was `marker_strong` until 2026-08-25, when Tim asked for it darker --
  # but that grey is also the top of the residual ramp, and the residual panel
  # was to stay exactly as it was. One value could not do both jobs.
  #
  # The cap is the end of the 95% interval, which is the number a reader takes
  # off this chart, so it is the darkest ink on it.
  marker_cap      = "#3D3D3D",
  plot_background = "#FFFFFF",
  grid_line       = "#D1D9E0",
  axis_line       = "#B0B8BF",
  # The A/E = 100% line. `ChartStyle` calls this the special line, and ships it
  # as the orange `#ffa742` with a medium blue commented out beside it. Tim
  # picked the blue on 2026-08-20 from the two rendered samples.
  #
  # LIGHTENED FROM THE SOURCES' `#1F77B4` ON 2026-08-25, in three passes against
  # rendered comparisons, along one line in CIELAB: L* 47.9, 58.3, 63.4, 68.7.
  #
  # The reason is not the reference line itself: at full strength it swallowed
  # the dark middle of any interval crossing it, so the one thing the line
  # exists to show -- whether a record's interval reaches 100% -- was the thing
  # it hid. Paler, the line can afford to be wider again, which is why the two
  # moved in opposite directions.
  reference_line  = "#6DB0D2",
  # The off-scale chevron. This is the orange `ChartStyle` ships for the special
  # line above, which Tim passed over for the blue -- so the value is the
  # prototype's and only its job is new. It has to carry on its own, being the
  # whole of what an off-scale record gets drawn, which is why it is the one
  # saturated colour on the chart apart from the 100% line.
  off_scale       = "#ffa742",
  # The residual panel's tint bands, `Background1To2` and `Background2To3` in
  # `AEChart.cs`. THE THREE SOURCES DISAGREE ON HOW STRONG THEY ARE: the C# says
  # 6.3% and 12.5%, `Test AE chart.svg` says 10% and 20%, and slide 6's
  # flattened fills work out at about 12.5% and 33%. The SVG is taken as the
  # most likely to be current, being an actual render rather than a constant or
  # a drawing.
  residual_band   = "#CF7D63",
  title_text      = "#1F2328",
  axis_title_text = "#3B424A",
  axis_scale_text = "#59636E"
)

#### Geometry ####
#
# All in pixels at 96 to the inch, which is how the prototype is written.
#
# The marker itself is drawn by `makeContent.aev_interval`, straight into
# `grid`, whose `lwd` is ALREADY in ninety-sixths of an inch and whose lengths
# can be given in inches. So those constants are used as they stand.
#
# The reference line and the theme go through ggplot2 instead, which wants its
# own unit, and `line_width()` is that conversion. It was MEASURED rather than
# read: rendering to an SVG device and reading the numbers back shows that a
# `linewidth` of 1 is `.pt` pixels, not the `96 / 25.4` that "linewidth is in
# millimetres" would imply. The difference is 33%.
pt_per_mm <- 72.27 / 25.4

# `grid` states a font size in big points, of which there are 72 to the inch
# against the design's 96 pixels.
pt_per_px <- 72 / 96

line_width <- function(pixels) pixels / pt_per_mm

# `CICircleRadius`, `CIThinLineWidth`, `CIThickLineWidth`, `CIBarHeight` and
# `CIBarWidth` in `AEChart.cs`, each checked against `Test AE chart.svg`, which
# is a rendering the prototype was drawn from: circles at `r="5"` with
# `stroke-width="0.75"`, arms at 2.5, and a cap `<line>` 18 long and 4 thick.
#
# Only the two weights Tim asked to lighten differ from those, and the note
# below says so. `cap_width_px` in particular is the design's own 18: it was
# briefly a fraction of a record's width, because an ordinary ggplot2 layer
# cannot be given a length in pixels, and the custom geom removed the need.
aev_geometry <- list(
  marker_radius_px = 5,
  # 0.75 IS ALSO THE CIRCLE'S STROKE, which is why it is left alone: the source
  # SVG draws the thin arms and the circle outline at one weight, and raising
  # this thickened the marker's outline along with them. The verticals Tim
  # called thin are the outer ones, which are `thick_width_px` below.
  thin_width_px    = 0.75,
  # WHAT A LINE WIDTH IS READ AGAINST IS ONE RECORD'S WIDTH, not the panel's.
  # That correction is the whole history of these two numbers.
  #
  # They were lightened to 0.75 and 1.75 on 2026-08-20, from the 2.5 the C# and
  # `Test AE chart.svg` both carry, on the reasoning that a pixel width reads
  # heavier the smaller the chart is drawn. True as far as it goes, and it was
  # measured against the panel -- but the panel is shared out among the records,
  # so the density that matters is pixels per record.
  #
  # `Test AE chart.svg` is 14 records across a 600 px panel: 43 px each, and a
  # 2.5 px arm is 5.8% of a record. Eight records at 8 by 7 inches gives 79 px
  # each, where the same arm would be 3.2% and the lightened 1.75 only 2.2%.
  # So the lightening and a roomier chart pulled the same way and compounded,
  # which is what Tim was seeing on 2026-08-25 when he called the verticals thin.
  #
  # RESTORED TO THE SOURCE VALUE, therefore, rather than lightened further.
  thick_width_px   = 2.5,
  # BACK TO THE SOURCE VALUE ON 2026-08-25, with the arm. Same reason: the
  # 2026-08-20 lightening was calibrated against the panel rather than against
  # one record's width. All three weights now match `Test AE chart.svg`.
  cap_height_px    = 4,
  cap_width_px     = 18,
  # `GridLineWidth` in `ChartStyle`, and `SpecialLineWidth` was 2.5 -- five
  # times it. SETTLED AT THREE TIMES ON 2026-08-25, after going 2.5, 1.5, 2 and
  # back to 1.5 across that day. The colour and the width were both cutting the
  # same problem -- the interval's dark middle being lost where it crossed -- so
  # the width came down while the blue was being lightened, went back up once it
  # was pale, and came down again when the arms were restored to their source
  # weight and no longer needed the line to hold its own against them.
  reference_width_px = 1.5,
  border_width_px  = 0.5,
  grid_width_px    = 0.5,
  # THE GAP ABOVE THE GROUP BAND, between the record names and the brackets
  # under them. Asked for on 2026-08-25: without it the names sit directly on
  # the bracket's end ticks and the two rows read as one block.
  #
  # It is a row of its own in the table rather than a margin on the strip text,
  # because the bracket is drawn to the edges of its cell -- a margin would move
  # the name and leave the verticals where they were.
  group_band_gap_px = 6,
  # The gap between groups. `DrawGroupClipRectanglesOpen` narrows each group by
  # at most `SmallMargin`, half at each side, so adjacent groups end up this far
  # apart.
  group_gap_px     = 6
)

# The off-scale chevron. NOTHING IN THE SOURCES HAS ONE: neither the pptx nor
# `AEChart.cs` marks a record it cannot draw, so every number here was settled
# by Tim from rendered mock-ups on 2026-08-25.
#
# It is deliberately blunt -- 32 wide against 7 high is an apex of 133 degrees,
# where a first attempt at 90 was hard to see. A chevron is not a scaled-down
# arrowhead and does not want to look like one.
aev_chevron <- list(
  width_px      = 32,
  height_px     = 7,
  line_width_px = 5,
  # How far the apex sits inside the panel edge, BEFORE the allowance for the
  # mitre that `makeContent.aev_chevron` adds. See the note there.
  inset_px      = 1
)

# The residual panel, from `DevianceResidualPlot` in `AEChart.cs`.
#
# ITS HEIGHT IS FIXED, NOT A SHARE. `MeasureVariableSize` returns a constant, so
# the strip is the same depth whatever the chart, and the A/E panel above takes
# whatever is left. Seven rows of 12 px, one per unit of residual from -3.5 to
# +3.5.
aev_residual <- list(
  limit             = 3.5,
  height_px         = 84,
  marker_radius_px  = 3,
  marker_width_px   = 1,
  # THE TOP OF THIS RAMP IS DELIBERATELY LOUD, raised from 2 on 2026-08-25 so
  # that a badly-fitting cell is unmissable rather than merely heavier. The
  # saturation point is unchanged: everything at or beyond 2.5 draws at 4, so a
  # residual of 9.69 and one of 2.6 still look alike and the printed number is
  # what separates them.
  drop_width_px     = c(0.5, 1, 4),
  # THE TOP OF THE RAMP IS AT 3, Tim's on 2026-08-25, up from the C#'s 2.5. It
  # now lands on the outer tint band's edge, so the residual at which a cell
  # stops getting louder is the one the bands stop marking.
  breakpoints       = c(0.5, 1.5, 3),
  band_alpha        = c(0.1, 0.2),
  # 9 px, which is the size of the record names under the chart. Tim's eye on
  # 2026-08-29, and deliberately NOT linked to `axis.text.x` in the theme: the
  # two happen to agree and are free to stop agreeing.
  text_px           = 9,
  # `TextGap` is a quarter of the marker radius in the C#, which put the value
  # almost against the zero line. WIDENED TO 3 ON 2026-08-25 at Tim's asking, so
  # that the numbers read as a row under the line rather than as part of it.
  text_gap_px       = 3,
  # How much smaller this panel's axis numbers and title are than the A/E
  # panel's. The title ratio is `DevYAxisTitleRatio`. The scale ratio is TIM'S,
  # asked for on 2026-08-21 and approved from a render -- the C# has 8 px
  # against the main axis's 9, which is barely a difference and not what he
  # meant by "much smaller".
  # THE TITLE IS NOT SHRUNK AT ALL, from 2026-08-25. It was 0.8, the C#'s own
  # DevYAxisTitleRatio, until Tim asked for the two panels' y-axis titles to
  # match. Only the numbers on the axis are made smaller now, which is what
  # "much smaller" meant in the first place.
  # THE TWO SCALES ARE THE SAME SIZE, which is what `ChartStyle` says: both get
  # `YAxisScaleFontSize`, and this ratio is what says so.
  #
  # It has moved twice and both moves are worth keeping. It was 0.7 from
  # 2026-08-21 -- "much smaller" -- which compounded with a smaller base to put
  # these numbers at 5.04 pt, under everything else on the chart. It went to
  # 10/9 on 2026-08-29 on the argument that this axis is always the same seven
  # values so nothing crowds it. Seven numbers in a fixed 84 px strip crowd
  # anyway, which is a thing only a render shows, so it came back to parity.
  scale_ratio       = 1,
  # `DevYAxisTitleRatio` in `AEChart.cs`. The design shrinks the residual
  # panel's title where it leaves the A/E panel's at full size, which is the
  # same distinction Tim asked for on 2026-08-25 before either of us had found
  # this constant.
  title_ratio       = 0.8
)

# THE THREE BLOCKS OF TEXT AT THE TOP AND BOTTOM, all from `ChartStyle`.
#
# `TitleFontSize` is 24 CSS px against the design's 12 px base text, so the
# title is TWICE the base size. ggplot2's own default is rel(1.2), which is
# where this chart sat until 2026-08-28 -- and the difference shows in the
# ratio a reader actually sees: the prototype's title is 2.7 times its axis
# text, where rel(1.2) made it 1.5 and the title read as timid.
#
# `DefaultMargin` is 10 px against the same 12 px base, so five sixths of the
# base size, and it is what slide 2 puts between the stacked blocks at the top.
# Expressed as ratios, all of it follows `base_size` like everything else.
aev_headings <- list(
  title_px    = 24,
  # THE DESIGN HAS NO SUBTITLE, so this one number is not `ChartStyle`'s. Tim
  # chose three quarters of the title on 2026-08-29, which puts it at 18 px --
  # above the panel titles at 16 and below the title at 24, which is the order
  # a reader needs.
  subtitle_px = 18,
  margin_ratio = 10 / 12
)

# NOTES ARE QUIETER THAN THE CHART, AND SET FURTHER APART THAN ORDINARY LINES.
#
# Measured before changing anything on 2026-08-28: the caption resolved to 7.2
# pt, exactly the size of the record names and the legend, so the notes competed
# with the chart's own labels instead of sitting under them. Worse, `lineheight`
# came through as theme_minimal's 0.9, which sets lines closer together than
# their own height -- fine for a wrapped sentence, wrong for a list where each
# line is a separate statement.
#
# So notes are smaller than anything else on the chart, and their lines are set
# well apart. All three numbers are ratios or multiples, so they follow the base
# size like everything else.
aev_notes <- list(
  # `XAxisScaleFontSize`, the smallest text the design has.
  text_px     = 9,
  # Of the notes' own size. 0.9 was cramped; 1.6 reads as separate notes.
  line_height = 1.6,
  # Between the chart and the first note, as a multiple of the base size.
  gap_ratio   = 1.0
)

# The key is drawn at the size the marker is drawn on the chart, so it has to
# be wide enough for a cap and tall enough for a circle with room either side.
aev_legend <- list(
  key_width_px  = 24,
  key_height_px = 16,
  text_ratio    = 0.8
)

#### Scale ####
#
# EVEN STEPS IN LOG A/E, NOT IN A/E. The axis is log A/E from -50% to +50% in
# steps of 10%, which comes out as 61, 67, 74, 82, 90, 100, 111, 122, 135, 149
# and 165 per cent. That is slide 3's table exactly, it is `DefaultMinLogAE` and
# `DefaultMaxLogAE` in `AEChart.cs`, and it is what the rendered samples show.
#
# Slides 6 and 7 disagree: they run 50 to 200 per cent in eighths of an octave.
# Tim ruled on 2026-08-20 that those slides were drawn before the range was
# settled and that this is the default.
#
# The axis runs from the first tick to the last with no expansion:
# `AEChart.cs` passes `ensureHalfUnitAtEnds: false` for the y axis, so the
# half-interval gap the layout rules ask for applies to the x axis only.
aev_log_ratio_limits <- c(-0.5, 0.5)
aev_log_step <- 0.1

# THE ENDS ARE SNAPPED, and it matters more than it looks. Repeated addition of
# the step lands a hair outside the range -- `c(-0.6, 0.6)` in steps of 0.2 ends
# at 0.6000000000000001 -- and ggplot2 censors an out-of-range break to NA,
# which silently loses the topmost grid line and its label.
aev_log_breaks <- function(log_range, log_step) {

  steps <- round(diff(log_range) / log_step)
  breaks <- log_range[1] + log_step * seq(0, steps)

  breaks[1] <- log_range[1]
  breaks[steps + 1L] <- log_range[2]
  breaks
}

aev_ratio_breaks <- exp(aev_log_breaks(aev_log_ratio_limits, aev_log_step))
# WHY THE RANGE IS BOUNDED AT ONE. Both panels share a single scale, and the
# only thing telling them apart is how far each reaches -- anything past 1 is
# the residual panel, anything inside it is A/E. A log A/E range of 1 is 37% to
# 272%, far wider than any of the three the design has ever used, so the bound
# costs nothing and keeps that test airtight.
#
# The step has to divide the range because the axis has no expansion: its ends
# ARE its outermost grid lines, and a step that fell short would leave the panel
# open at the top. Dividing both ends also puts a break on 100% for free, which
# matters rather more -- that is the line everything is read against.
aev_check_log_scale <- function(log_range, log_step) {

  if (!is.numeric(log_range) || length(log_range) != 2L || anyNA(log_range) ||
        log_range[1] >= log_range[2]) {
    stop("`log_range` must be two increasing numbers.", call. = FALSE)
  }

  if (max(abs(log_range)) >= 1) {
    stop("`log_range` must lie within -1 and 1, which is an A/E of 37% to ",
         "272%. Beyond that the A/E panel could not be told from the residual ",
         "panel, which is how one scale serves both.", call. = FALSE)
  }

  if (!is.numeric(log_step) || length(log_step) != 1L || is.na(log_step) ||
        log_step <= 0) {
    stop("`log_step` must be a single positive number.", call. = FALSE)
  }

  divides <- function(value) abs(value / log_step - round(value / log_step)) < 1e-8

  if (!all(vapply(log_range, divides, logical(1)))) {

    # The nearest range that would work, rounded outwards. Saying so matters:
    # the most natural request of all -- twenty per cent steps on the default
    # range -- is refused, because 0.2 does not divide 0.5.
    widened <- c(
      floor(log_range[1] / log_step) * log_step,
      ceiling(log_range[2] / log_step) * log_step
    )

    stop("`log_step` must divide both ends of `log_range`, so that the axis ",
         "ends on a labelled tick and 100% is one of them. With `log_step` of ",
         format(log_step), ", a `log_range` of c(", format(widened[1]), ", ",
         format(widened[2]), ") would work.", call. = FALSE)
  }

  invisible(NULL)
}

# 1.96, to the precision `[[.aev` uses for `log_A_over_E_95pc`.
aev_z_975 <- 1.9599639845400540

# NOTES ARE A LIST OF STRINGS AND NOTHING MORE. Tim's rule, 2026-08-28.
#
# The chart is not coupled to `aev()` -- an `aev` is three vectors and anybody
# may build one -- so the notes must not be coupled to it either. Nothing here
# knows what an overdispersion is, or a Z, or a test mortality. It takes lines
# of text and puts them under the chart, which is the one shape that ports:
# a list of strings is expressible in R, in Python and in C# alike, and needs no
# type of its own in any of them.
#
# One string per line, in the order given. Empty lines are kept, so a caller can
# space a block of notes if they want to.
aev_notes_text <- function(notes) {

  if (length(notes) == 0L) return(NULL)

  if (!is.character(notes)) {
    stop("The argument `notes` must be a character vector, one line per note.",
         call. = FALSE)
  }

  if (anyNA(notes)) {
    stop("The argument `notes` cannot contain missing values.", call. = FALSE)
  }

  paste(notes, collapse = "\n")
}

# The two panels. These strings are the y-axis titles: `facet_grid(switch = "y")`
# puts the strip on the left of the panel it belongs to, which is exactly where
# slide 2's layout stack puts a y-axis title, and is the only way to give two
# panels two different ones inside a single plot.
#
# THE QUALIFIER GOES ON ITS OWN LINE, so the name reads first and the gloss
# second, and neither title takes more width than it needs.
#
# "(+/-2 ~ 95%)" IS THE ONLY THING ON THE CHART THAT CONNECTS THE TWO PANELS.
# The upper one is labelled 95% confidence and the lower one is the exact test
# that interval approximates, and without this nothing says so. Tim put it in
# the title rather than the legend, which is better: it sits against the panel
# it explains and survives a chart drawn without a legend.
#
# The non-breaking space is `Title_AEOnLogScale` in `AEChart.cs`, kept so the
# qualifier cannot come apart. Written as escapes so the file stays ASCII.
aev_panels <- c(
  "A/E\n(log\u00A0scale)",
  "Deviance residual\n(\u00B12 ~ 95%)"
)

# BOTH PANELS ARE PLOTTED IN LOG A/E AND IN RESIDUALS RESPECTIVELY, on one
# untransformed scale. A transformed scale is not available here: it would take
# the logarithm of the residuals too.
#
# So one scale has to label two panels differently, and it dispatches on the
# range, which is safe because both ranges are constants of this chart and are
# an order of magnitude apart. Testing against the residual limit rather than
# for equality matters: the range reaching these functions may carry the
# scale's expansion, and an equality test silently fell through to the A/E
# labels when it did.
aev_panel_limits <- function(range, log_range = aev_log_ratio_limits) {
  if (max(abs(range)) > 1) c(-aev_residual$limit, aev_residual$limit) else log_range
}

aev_panel_breaks <- function(limits, log_breaks = log(aev_ratio_breaks)) {
  if (max(abs(limits)) > 1) seq(-3, 3) else log_breaks
}

# The same test as the other two rather than `breaks %in% -3:3`, which could be
# satisfied by an A/E axis whose own breaks happened to be whole numbers.
# `na.rm`, because ggplot2 hands this the breaks with any it censored replaced
# by NA, and one NA would otherwise decide nothing at all.
aev_panel_labels <- function(breaks) {
  # Plotmath for the residual panel, plain text for the A/E panel: percentages
  # have no sign to draw and their `%` would not parse.
  if (max(abs(breaks), na.rm = TRUE) > 1) return(aev_signed_math(breaks, 0))
  paste0(round(exp(breaks) * 100), "%")
}

# Slide 1: "Always use minus sign for displaying numbers, never hyphen dash."
# Zero takes no sign at all, as the rendered samples have it.
# THE SAME NUMBERS AS PLOTMATH, WHICH IS THE ONLY PORTABLE WAY TO DRAW THEM.
#
# `aev_signed()` above uses U+2212, the true minus, because that is what the
# design asks for and it is right in a string. It cannot be DRAWN portably: a
# plain string is converted to the device's encoding first, U+2212 is not in
# latin1, and `R CMD check` runs in a latin1 locale. The chart errored out with
#
#     conversion failure in 'mbcsToSbcs': for <U+2212>
#
# Plotmath never converts. It selects the SYMBOL font and draws a glyph from it,
# where byte 0x2D is a true minus rather than the hyphen the same byte means in
# Helvetica. Measured from an uncompressed PDF on 2026-08-29:
#
#     "-3"             /F2 (-3) Tj                  one font, a hyphen
#     expression(-3)   /F6 (-) Tj  /F2 (3) Tj       F6 is /BaseFont /Symbol
#
# The plus comes from Symbol too, at the same position and the same width, so
# the two signs match and a column of signed residuals stays aligned. That was
# Tim's question and is the reason no separate handling is needed for the plus.
#
# THE DIGITS ARE QUOTED AND THE SIGN IS NOT, which is the whole trick.
# `parse(text = "-2.60")` renders "2.6" -- plotmath treats the digits as a
# number and drops the trailing zero, silently reformatting every residual on
# the chart. `parse(text = '-"2.60"')` renders the sign as an operator, in
# Symbol, and the digits as a string, in the text font, exactly as formatted.
aev_signed_math <- function(values, digits) {

  formatted <- formatC(abs(values), format = "f", digits = digits)
  sign <- ifelse(values > 0, "+", ifelse(values < 0, "-", ""))

  parse(text = paste0(sign, "\"", formatted, "\""))
}

#### Status ####

# WHY THIS IS NOT A COLUMN OF `as.data.frame.aev()`. A chart has to tell four
# cases apart, and only three of them are properties of the triple:
#
# * `ok`       -- has a position on both panels.
# * `empty`    -- A, E and V all zero. A true statement about nobody.
# * `missing`  -- any of A, E or V missing. An absence of information.
# * `off_scale`-- a real answer with nowhere to put it, either because A is zero
#                 and log(A/E) is therefore infinite, or because A/E is outside
#                 the axis.
#
# The fourth depends on the axis, which is a property of the chart. That is the
# whole reason the column was kept out of the frame and the work put here.
#
# `empty` and `missing` are indistinguishable in every calculated property --
# both read NaN -- so this reads the raw triple and nothing else.
# LIMITS ARE IN LOG A/E, not in A/E. They used to be a ratio and were logged
# again inside, which put this test and the chevron's on opposite sides of an
# `exp()`/`log()` round trip -- and that round trip is not exact: 0.3, 0.35 and
# 0.45 all come back a unit in the last place short, so a record sitting exactly
# on such a limit could be off scale here and unmarked there, and be drawn
# nowhere at all. Found by review on 2026-08-26.
aev_status <- function(x, limits = aev_log_ratio_limits) {

  u <- unclass(x)
  A <- u[["A"]]
  E <- u[["E"]]
  V <- u[["V"]]

  ratio <- A / E

  # ONE STANDARD DEVIATION, NOT THE POINT AND NOT THE 95% ARMS. Tim's rule,
  # 2026-08-25, having seen a chart drawn the other way: a record belongs on the
  # chart if the marker and its inner interval reach it, and a record that does
  # not is replaced by a chevron even where its outer arms would have shown.
  #
  # Drawing the outer arms of a marker whose centre is far off the top gives a
  # bar with no dot and no cap, which says less than the chevron does and takes
  # more ink to say it.
  centre <- log(ratio)
  reach <- sqrt(V) / E
  top <- limits[2]
  bottom <- limits[1]

  # `%in% TRUE` rather than a bare comparison: a NaN operand gives NA, and NA
  # cannot be used as a subscript in an assignment.
  is_empty <- (A == 0 & E == 0 & V == 0) %in% TRUE
  is_missing <- is.na(A) | is.na(E) | is.na(V)
  is_off_scale <- !((centre + reach >= bottom & centre - reach <= top) %in% TRUE)

  status <- rep("ok", length(A))
  status[is_off_scale] <- "off_scale"
  status[is_empty] <- "empty"
  status[is_missing] <- "missing"

  factor(status, levels = c("ok", "empty", "missing", "off_scale"))
}

# One row per record, with everything the panels need. Records with no labels
# are numbered, so the x axis always says which record is which.
aev_plot_frame <- function(x, limits = aev_log_ratio_limits) {

  frame <- as.data.frame(x)
  frame$position <- seq_len(nrow(frame))
  frame$status <- aev_status(x, limits)
  frame$log_ratio <- log(frame$A_over_E)

  # WHAT THE SCALE IS ALLOWED TO SEE. See `GeomAevInterval`: a drawable record
  # may sit well outside the axis, and training the scale on it turns the A/E
  # panel into a residual one.
  frame$clamped <- pmin(pmax(frame$log_ratio, limits[1]), limits[2])

  unlabelled <- is.na(frame$name)
  frame$name[unlabelled] <- as.character(frame$position[unlabelled])

  frame
}

#### The marker ####

# ONE GROB FOR THE WHOLE MARKER, because the design mixes two kinds of length
# and an ordinary layer can only speak one of them.
#
# Where the interval ends is a fact about the data: the arms run to one and to
# 1.96 standard deviations of log A/E. How big the circle is, how wide the cap
# is, and where the thin line has to stop are facts about the picture, fixed in
# pixels whatever the data says. A `geom_linerange` can only be given data
# units, which is why the cap width had to be guessed as a fraction of a record
# and why the thin line could not be stopped at the circle at all.
#
# `grid` resolves both at draw time, when the panel size is finally known. So
# this hands it the data-space positions as `npc` and the pixel sizes as
# inches, and lets it do the arithmetic. Two consequences worth knowing:
#
# * `grid`'s `lwd` is already in ninety-sixths of an inch, i.e. in the pixels
#   the design is written in, so no conversion happens here at all.
# * The circle is a real circle of a stated radius, not a plotting character
#   whose diameter has to be worked back out of a font size.
#
# EVERY LINE IS DRAWN WITH BUTT ENDS. `grid` rounds line ends by default, and
# `AEChart.cs` draws these as rectangles, so the default gave the arms and caps
# rounded corners. It also lengthened them: a round end adds half a line width
# at each end, so an 18 px cap 3 px thick was covering 21 px.

geom_aev_interval <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::layer(
    geom = GeomAevInterval,
    mapping = mapping,
    data = data,
    stat = "identity",
    position = "identity",
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(...)
  )
}

GeomAevInterval <- ggplot2::ggproto(
  "GeomAevInterval",
  ggplot2::Geom,

  # `sigma` is the standard deviation of log A/E. It is not a position, so no
  # scale touches it, which is what we want: `y` gets the log transform applied
  # to it and `sigma` is already a length in that transformed space.
  # `y` AND `centre` ARE BOTH THE SAME NUMBER FOR ALMOST EVERY RECORD, and the
  # difference is the whole point. `y` is log A/E CLAMPED into the axis, and is
  # what the scale is trained on; `centre` is where the marker is actually
  # drawn. Nothing here reads `y`.
  #
  # WITHOUT THE SPLIT THE A/E PANEL SILENTLY BECOMES A RESIDUAL PANEL. A record
  # is drawable when its inner interval reaches the axis, which its centre need
  # not do -- so an ordinary small cell at A/E = 300% with a wide interval is
  # drawn, and trains the panel's range out to 1.1. The two panels are told
  # apart by how far their ranges reach, so the A/E panel was then read as the
  # residual one and drawn on the deviance axis, labels and all, with no error.
  # Found by review on 2026-08-26; `test-aev-scale.R` holds the line now.
  required_aes = c("x", "y", "centre", "sigma"),
  default_aes = ggplot2::aes(),
  draw_key = ggplot2::draw_key_blank,

  draw_panel = function(data, panel_params, coord) {

    if (nrow(data) == 0L) return(grid::nullGrob())

    # Each y of interest through the coordinate system separately. Passing the
    # whole thing through once would not do: `coord$transform` only knows how
    # to place a column called `y`.
    place <- function(values) {
      coord$transform(data.frame(x = data$x, y = values), panel_params)$y
    }

    positions <- coord$transform(data, panel_params)$x

    grid::gTree(
      x           = positions,
      y           = place(data$centre),
      inner_lower = place(data$centre - data$sigma),
      inner_upper = place(data$centre + data$sigma),
      outer_lower = place(data$centre - aev_z_975 * data$sigma),
      outer_upper = place(data$centre + aev_z_975 * data$sigma),
      cl = "aev_interval"
    )
  }
)

# Pixels to a length grid understands. The design is written at 96 per inch.
inches <- function(pixels) grid::unit(pixels / 96, "inches")

#' @importFrom grid makeContent
# HOW BIG A MARK IS, RELATIVE TO THE TYPE AROUND IT.
#
# Every pixel length on a panel is a multiple of the chart's base text size, so
# that a chart drawn larger is resized by asking for larger type rather than by
# a second knob that has to be kept in step with the first. That is ggplot2's
# own convention -- theme_grey() derives base_line_size from base_size -- and at
# the reference size below every constant reproduces its old pixel exactly, so
# adopting it changed no output at all.
#
# The reference is 9 because that is theme_aev()'s default, and the size every
# constant in this file was drawn at.
aev_reference_text_size <- 9

# WHERE THE SIZE COMES FROM, and why it is not an argument to autoplot().
#
# A base_size argument would be the obvious way, and would be wrong: the marks
# would be fixed when the plot was made, and a later + theme(...) would move the
# type without them. grid cannot see the theme at draw time -- but the FACET is
# handed the finished theme while the table is assembled, which is after every
# + has been applied and before any mark has resolved its lengths. So the facet
# reads the size and stamps it on, and the marks follow what the theme says.
aev_mark_scale <- function(theme) {

  if (is.null(theme)) return(1)

  size <- ggplot2::calc_element("text", theme)$size

  if (is.null(size) || !is.finite(size) || size <= 0) return(1)

  size / aev_reference_text_size
}

# Named pixel lengths multiplied by the scale, everything else left alone. The
# fields are named at each call rather than scaled wholesale, because these
# lists hold ratios and residual sizes as well as lengths.
aev_scaled <- function(constants, scale, fields) {

  if (scale == 1) return(constants)

  constants[fields] <- lapply(constants[fields], function(value) value * scale)
  constants
}

# The marks that carry a scale. Each is a gTree that resolves its own lengths at
# draw time, which is what makes stamping a field on it work at all.
aev_scalable_marks <- c("aev_interval", "aev_residual_marker", "aev_chevron")

aev_stamp_scale <- function(grob, scale) {

  if (inherits(grob, aev_scalable_marks)) {
    grob$scale <- scale
    return(grob)
  }

  # Children are edited in place: rebuilding the list with gList() would drop
  # the names a gTree matches against its childrenOrder.
  for (index in seq_along(grob$children)) {
    grob$children[[index]] <- aev_stamp_scale(grob$children[[index]], scale)
  }

  for (index in seq_along(grob$grobs)) {
    grob$grobs[[index]] <- aev_stamp_scale(grob$grobs[[index]], scale)
  }

  grob
}

#' @exportS3Method grid::makeContent
makeContent.aev_interval <- function(x) {

  geometry <- aev_scaled(
    aev_geometry, x$scale %||% 1,
    c("marker_radius_px", "thin_width_px", "thick_width_px",
      "cap_height_px", "cap_width_px")
  )
  palette <- aev_palette

  # THE ONE THING THAT CANNOT BE DECIDED ANY EARLIER: whether the interval is
  # longer than the circle is wide. One is a data length and the other is a
  # pixel length, and they only become comparable once there is a viewport to
  # measure against. Here there is one.
  radius <- grid::convertHeight(inches(geometry$marker_radius_px), "npc", valueOnly = TRUE)
  half_cap <- grid::convertWidth(inches(geometry$cap_width_px / 2), "npc", valueOnly = TRUE)

  # THE THIN LINE STOPS AT THE CIRCLE and is left out altogether when the
  # circle already covers it. The thick arms are not treated this way: the
  # design says they are drawn even where they intrude.
  showing <- abs(x$inner_upper - x$y) > radius

  thin <- if (any(showing)) {
    grid::segmentsGrob(
      x0 = grid::unit(rep(x$x[showing], 2L), "npc"),
      x1 = grid::unit(rep(x$x[showing], 2L), "npc"),
      y0 = grid::unit(c(x$y[showing] + radius, x$y[showing] - radius), "npc"),
      y1 = grid::unit(c(x$inner_upper[showing], x$inner_lower[showing]), "npc"),
      gp = grid::gpar(col = palette$marker_weak, lwd = geometry$thin_width_px, lineend = "butt")
    )
  } else {
    grid::nullGrob()
  }

  thick <- grid::segmentsGrob(
    x0 = grid::unit(rep(x$x, 2L), "npc"),
    x1 = grid::unit(rep(x$x, 2L), "npc"),
    y0 = grid::unit(c(x$inner_upper, x$inner_lower), "npc"),
    y1 = grid::unit(c(x$outer_upper, x$outer_lower), "npc"),
    gp = grid::gpar(col = palette$marker_medium, lwd = geometry$thick_width_px, lineend = "butt")
  )

  caps <- grid::segmentsGrob(
    x0 = grid::unit(rep(x$x, 2L) - half_cap, "npc"),
    x1 = grid::unit(rep(x$x, 2L) + half_cap, "npc"),
    y0 = grid::unit(c(x$outer_upper, x$outer_lower), "npc"),
    y1 = grid::unit(c(x$outer_upper, x$outer_lower), "npc"),
    gp = grid::gpar(col = palette$marker_cap, lwd = geometry$cap_height_px, lineend = "butt")
  )

  circle <- grid::circleGrob(
    x = grid::unit(x$x, "npc"),
    y = grid::unit(x$y, "npc"),
    r = inches(geometry$marker_radius_px),
    gp = grid::gpar(
      col = palette$marker_weak,
      # Transparent, so the reference line and the grid show through.
      fill = NA,
      lwd = geometry$thin_width_px
    )
  )

  # Slide 7's order: circle, thin, thick, caps.
  grid::setChildren(x, grid::gList(circle, thin, thick, caps))
}

#### Records that are off the scale entirely ####

# WHICH WAY AN OFF-SCALE RECORD WENT. Whether it is off the scale at all is
# `aev_status()`'s question and is asked the same way there -- on the inner
# interval, A/E give or take one standard deviation -- so the two must and do
# agree. This only sorts them into the two edges.
#
# TWO CASES FALL OUT WITHOUT BEING WRITTEN. A record with no deaths has a
# centre of -Inf and a finite reach, so -Inf plus the reach is still -Inf and it
# tests as below the bottom like any other -- which is Tim's ruling that A = 0
# is treated the same as anything else off the scale, obtained for nothing.
# `empty` and `missing` both read NaN, every comparison gives NA, and `%in% TRUE`
# turns that into FALSE, so neither is marked. That is right rather than merely
# convenient: a chevron says which way a record went, and neither "nobody" nor
# "not known" went anywhere.
aev_off_scale <- function(frame, log_range) {

  # WHETHER A RECORD IS OFF SCALE IS ALREADY DECIDED, by `aev_status()`. This
  # only sorts those it flagged into the two edges, and deliberately does not
  # ask the question again: asking it twice is what let the two answers drift
  # apart, and a record they disagreed about was drawn nowhere at all.
  off <- frame$status == "off_scale"
  centre <- frame$log_ratio

  above <- off & ((centre > log_range[2]) %in% TRUE)
  below <- off & ((centre < log_range[1]) %in% TRUE)

  # The whole row is kept, not just the position: the layer needs `panel` and
  # `group` or the facet draws it in every column.
  chevrons <- rbind(frame[above, , drop = FALSE], frame[below, , drop = FALSE])
  chevrons$top <- c(rep(TRUE, sum(above)), rep(FALSE, sum(below)))

  chevrons
}

geom_aev_chevron <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::layer(
    geom = GeomAevChevron,
    mapping = mapping,
    data = data,
    stat = "identity",
    position = "identity",
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(...)
  )
}

GeomAevChevron <- ggplot2::ggproto(
  "GeomAevChevron",
  ggplot2::Geom,

  # `top` says which edge to sit against. Like `sigma` on the marker it is not
  # a position, so no scale touches it.
  required_aes = c("x", "top"),
  default_aes = ggplot2::aes(),
  draw_key = ggplot2::draw_key_blank,

  draw_panel = function(data, panel_params, coord) {

    if (nrow(data) == 0L) return(grid::nullGrob())

    # A CHEVRON HAS NO Y. It sits against the edge of the panel whatever the
    # record's value is, and the value is off the scale in any case. But
    # `coord$transform` places a column called `y` and refuses without one, so
    # it is given a zero -- inside every allowed range, so nothing is censored
    # -- and the transformed y is thrown away.
    positions <- coord$transform(data.frame(x = data$x, y = 0), panel_params)$x

    grid::gTree(x = positions, top = data$top, cl = "aev_chevron")
  }
)

#' @exportS3Method grid::makeContent
makeContent.aev_chevron <- function(x) {

  chevron <- aev_scaled(
    aev_chevron, x$scale %||% 1,
    c("width_px", "height_px", "inset_px", "line_width_px")
  )

  half <- grid::convertWidth(
    inches(chevron$width_px), "npc", valueOnly = TRUE
  ) / 2
  height <- grid::convertHeight(
    inches(chevron$height_px), "npc", valueOnly = TRUE
  )

  # THE INSET HAS TO CLEAR THE MITRE. A mitred join puts ink half a line width
  # beyond the apex point, so an inset measured to the point alone clips the tip
  # off anything thicker than twice it -- which the 5 px chevron is.
  inset <- grid::convertHeight(
    inches(chevron$inset_px + chevron$line_width_px / 2),
    "npc", valueOnly = TRUE
  )

  apex <- ifelse(x$top, 1 - inset, inset)
  base <- ifelse(x$top, 1 - inset - height, inset + height)

  # One polyline of three points rather than two segments: a mitred join keeps
  # the apex sharp, where two butt-ended segments would leave a notch in it.
  chevrons <- grid::polylineGrob(
    x = grid::unit(as.vector(rbind(x$x - half, x$x, x$x + half)), "npc"),
    y = grid::unit(as.vector(rbind(base, apex, base)), "npc"),
    id = rep(seq_along(x$x), each = 3),
    gp = grid::gpar(
      col = aev_palette$off_scale,
      lwd = chevron$line_width_px,
      lineend = "butt",
      linejoin = "mitre"
    )
  )

  grid::setChildren(x, grid::gList(chevrons))
}

#### The residual panel ####

# HOW LOUD A RESIDUAL LOOKS IS A FUNCTION OF ITS SIZE, which is the panel's one
# real idea: the drop line thickens and darkens as the residual grows, so a bad
# cell is visible before any number is read. `AEChart.cs` sets the breakpoints
# at 0.5, 1.5 and 2.5 and interpolates between them, in CIELAB for the colour
# because that is where interpolation looks even.
#
# Above 2.5 nothing more happens. A residual of 6 draws like a residual of 3,
# and the printed number is what separates them.
aev_residual_ramp <- function(z, low, mid, high) {

  magnitude <- abs(z)
  breaks <- aev_residual$breakpoints

  first <- pmin(1, pmax(0, (magnitude - breaks[1]) / (breaks[2] - breaks[1])))
  second <- pmin(1, pmax(0, (magnitude - breaks[2]) / (breaks[3] - breaks[2])))

  ifelse(magnitude <= breaks[2],
         low + first * (mid - low),
         mid + second * (high - mid))
}

aev_residual_width <- function(z) {
  widths <- aev_residual$drop_width_px
  aev_residual_ramp(z, widths[1], widths[2], widths[3])
}

# The same ramp on colour, done in CIELAB to match `LerpLab`.
aev_residual_colour <- function(z) {

  position <- aev_residual_ramp(z, 0, 0.5, 1)

  ramp <- grDevices::colorRamp(
    c(aev_palette$marker_weak, aev_palette$marker_medium, aev_palette$marker_strong),
    space = "Lab"
  )

  rgb <- ramp(position)
  interpolated <- grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)

  # THE THREE ANCHORS ARE THE PALETTE'S OWN GREYS, not what comes back from a
  # round trip into CIELAB and out again -- that loses a unit in each channel,
  # so a residual of 3 came out `#4C4C4C` where the palette says `#4D4D4D`.
  # Only the colours between them are computed.
  anchors <- c(aev_palette$marker_weak, aev_palette$marker_medium, aev_palette$marker_strong)
  at <- match(position, c(0, 0.5, 1))

  ifelse(is.na(at), interpolated, anchors[at])
}

geom_aev_residual <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::layer(
    geom = GeomAevResidual,
    mapping = mapping,
    data = data,
    stat = "identity",
    position = "identity",
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(...)
  )
}

GeomAevResidual <- ggplot2::ggproto(
  "GeomAevResidual",
  ggplot2::Geom,

  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(),
  draw_key = ggplot2::draw_key_blank,

  draw_panel = function(data, panel_params, coord) {

    if (nrow(data) == 0L) return(grid::nullGrob())

    place <- function(values) {
      coord$transform(data.frame(x = data$x, y = values), panel_params)$y
    }

    grid::gTree(
      x = coord$transform(data, panel_params)$x,
      y = place(data$y),
      zero = place(rep(0, nrow(data))),
      # The residual itself, kept because the colour, the line width and the
      # printed label all depend on it and none of them survives the transform.
      z = data$y,
      cl = "aev_residual_marker"
    )
  }
)

#' @exportS3Method grid::makeContent
makeContent.aev_residual_marker <- function(x) {

  scale <- x$scale %||% 1
  residual <- aev_scaled(
    aev_residual, scale,
    c("marker_radius_px", "marker_width_px", "text_px", "text_gap_px")
  )

  radius <- grid::convertHeight(inches(residual$marker_radius_px), "npc", valueOnly = TRUE)
  gap <- grid::convertHeight(inches(residual$text_gap_px), "npc", valueOnly = TRUE)

  colours <- aev_residual_colour(x$z)

  # The ramp is a function of the residual and knows nothing of the chart, so
  # its widths are scaled here rather than inside it.
  widths <- aev_residual_width(x$z) * scale
  above <- x$z >= 0

  # SAME RULE AS THE THIN LINE ON THE MARKER: the drop line runs from the zero
  # line to the edge of the circle, and is left out where the circle already
  # covers the distance.
  reach <- abs(x$y - x$zero)
  showing <- reach > radius

  drop <- if (any(showing)) {
    grid::segmentsGrob(
      x0 = grid::unit(x$x[showing], "npc"),
      x1 = grid::unit(x$x[showing], "npc"),
      y0 = grid::unit(x$zero[showing], "npc"),
      y1 = grid::unit(
        ifelse(above[showing], x$y[showing] - radius, x$y[showing] + radius), "npc"
      ),
      gp = grid::gpar(col = colours[showing], lwd = widths[showing], lineend = "butt")
    )
  } else {
    grid::nullGrob()
  }

  circle <- grid::circleGrob(
    x = grid::unit(x$x, "npc"),
    y = grid::unit(x$y, "npc"),
    r = inches(residual$marker_radius_px),
    gp = grid::gpar(
      col = aev_palette$marker_medium,
      fill = NA,
      lwd = residual$marker_width_px
    )
  )

  # THE NUMBER GOES ON THE FAR SIDE OF THE ZERO LINE FROM ITS MARKER. That is
  # what makes the slide's claim that the values are "guaranteed to fit" true:
  # the marker can be anywhere in its half of the panel, and the other half is
  # always empty at that x. The extra offset only bites for a residual so small
  # that its circle straddles the line.
  crowding <- pmax(0, radius - reach)
  labels <- grid::textGrob(
    label = aev_signed_math(x$z, 2),
    x = grid::unit(x$x, "npc"),
    y = grid::unit(ifelse(above, x$zero - gap - crowding, x$zero + gap + crowding), "npc"),
    hjust = 0.5,
    vjust = ifelse(above, 1, 0),
    # ONE COLOUR FOR EVERY VALUE, whatever its residual. Tim reversed the
    # earlier ruling on 2026-08-25: the drop line already carries the size
    # through its weight and colour, and a number that faded with it was hard
    # to read in exactly the cells that were fine.
    gp = grid::gpar(
      col = aev_palette$marker_strong,
      fontsize = residual$text_px * pt_per_px
    )
  )

  grid::setChildren(x, grid::gList(drop, circle, labels))
}

#### Two panels, one x axis, between them ####

# WHAT THIS IS FOR, and why it is not a theme.
#
# Slide 2's layout stack runs plot area, x scale, group names, residual panel:
# the category labels sit BETWEEN the two charts, not under both. Both rendered
# samples do the same. `facet_grid` puts the shared x axis at the bottom, and a
# theme cannot say otherwise because an axis belongs to the layout rather than
# to the styling.
#
# Nor can a theme make one panel's axis numbers smaller than another's, which
# is the second thing wanted here: `axis.text.y` is one element for the whole
# plot.
#
# So both are done by subclassing the facet and editing the table it builds.
# Custom facets are an extension point rather than a trick -- `lemon` does the
# same to repeat axes -- and the important consequence is that the result is
# still an ordinary `ggplot` that a user can extend with `+`. Reaching for
# `ggplotGrob()` after the fact would give back a `gtable`, which they could not.
#
# HOW BRITTLE THIS IS, honestly. Two ggplot2 internals are relied on: that
# `axes = "all_x"` draws an x axis under every panel, and that the assembled
# table names its cells `panel`, `axis-b`, `axis-l` and `strip-l`. Those names
# have been stable for years and most of the ggplot2 ecosystem depends on them.
# The code below leans on them as little as it can -- it matches only those
# four prefixes and finds everything else by ROW POSITION, so a renumbering
# cannot mislead it. Every effect is asserted in `test-aev-facet.R`, so a
# ggplot2 change fails the suite rather than quietly producing a wrong chart.

# Numbering the runs. A record starts a new run when its group name differs from
# the one before it, and two records with no name at all are in the same run.
aev_runs <- function(names) {

  count <- length(names)
  if (count == 0L) return(factor(character(0)))

  previous <- names[-count]
  current <- names[-1]

  # NA HAS TO BE COMPARED BY HAND, because `NA != NA` is NA and cannot start a
  # run. Two unnamed records continue one; a name after an absence starts one.
  changed <- c(
    TRUE,
    is.na(previous) != is.na(current) |
      (!is.na(previous) & !is.na(current) & previous != current)
  )

  factor(cumsum(changed), levels = seq_len(sum(changed)))
}

facet_aev <- function(groups, label_lines = 1L, run_names = character(0)) {

  # GROUPS ARE A SECOND FACET DIMENSION, and that is what draws the gaps. A
  # column per group means a gap wherever the group changes and nowhere else,
  # which is the prototype's rule exactly -- and because a column spans the
  # whole grid, the gap runs through the residual panel too, which is the other
  # half of the rule.
  #
  # `space = "free_x"` makes each column as wide as the records it holds, so a
  # record is the same width whichever group it is in.
  base <- if (groups) {
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$panel),
      cols = ggplot2::vars(.data$run),
      # The run is a number; the bracket under it says the group's name. Two
      # runs of the same name are two brackets, both saying it.
      #
      # NAMED, through `labeller()`, because a bare function is handed the ROW
      # variable as well as the column one and would try to read the panel names
      # as run numbers.
      labeller = ggplot2::labeller(
        run = function(values) run_names[as.integer(as.character(values))]
      ),
      scales = "free",
      space = "free_x",
      switch = "both",
      axes = "all_x"
    )
  } else {
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$panel),
      scales = "free_y",
      switch = "y",
      # An x axis under EVERY panel. The one under the lower panel is removed
      # below; the one under the upper panel is the one the design asks for.
      axes = "all_x"
    )
  }

  ggplot2::ggproto(
    NULL, FacetAev,
    shrink = base$shrink, params = base$params, label_lines = label_lines
  )
}

# ggplot2 CALLS draw_panels POSITIONALLY -- names() on the dots comes back
# empty, checked on 4.0.3 -- so an argument cannot be picked out by name. The
# names are taken from the parent's own formals rather than written out here,
# which is the same reasoning as passing ... through to it: if a future ggplot2
# adds or reorders an argument, this follows it instead of breaking.
# WHEN A LOOKUP INTO ggplot2'S OWN LAYOUT COMES BACK EMPTY, STOP.
#
# The chart finds its way around the table ggplot2 builds by cell name --
# `panel`, `axis-b`, `strip-l`, `strip-b` -- and reads one argument of
# `FacetGrid$draw_panels` by position. Those are ggplot2's, not ours. They have
# been stable for years and much of the extension ecosystem leans on them, but
# a release could change any of them.
#
# THE POINT OF ERRORING IS THAT THE ALTERNATIVE IS WORSE. Every one of these
# lookups used to fail QUIETLY: a missed `axis-b` left the record names printed
# twice, a missed `strip-l` left the residual axis at full size, a missed
# `theme` left every mark at reference size while the type scaled. In each case
# the chart still drew, still looked plausible, and was wrong. An error naming
# what could not be found is recoverable; a confidently wrong chart is not.
#
# Ruled by Tim on 2026-08-26: loud for the structural lookups.
aev_require_layout <- function(found, what) {

  # `is.atomic` guards the NA test: a row number that came back NA is a miss,
  # but a THEME is a list holding plenty of legitimate NAs -- `fill = NA` on the
  # legend box, for one -- and testing that would refuse every chart.
  missing <- length(found) == 0L || (is.atomic(found) && anyNA(found))
  if (!missing) return(invisible(found))

  stop(
    "logmu could not find ", what, " in the plot that ggplot2 built.\n",
    "This usually means the installed ggplot2 has changed something the ",
    "logmu chart depends on.\n",
    "Installed ggplot2: ", as.character(utils::packageVersion("ggplot2")), ".",
    call. = FALSE
  )
}

aev_facet_arguments <- function(...) {

  given <- list(...)

  method <- ggplot2::FacetGrid$draw_panels
  inner <- environment(method)$f
  expected <- setdiff(names(formals(inner %||% method)), "self")

  supplied <- names(given)
  if (is.null(supplied)) supplied <- rep("", length(given))

  blank <- !nzchar(supplied)
  supplied[blank] <- expected[seq_along(given)][blank]

  stats::setNames(given, supplied)
}

FacetAev <- ggplot2::ggproto(
  "FacetAev",
  ggplot2::FacetGrid,

  # `...` straight through to the parent, deliberately: the argument list of
  # `draw_panels` has changed between ggplot2 versions and naming them here
  # would tie this to one of them.
  draw_panels = function(self, ...) {

    table <- ggplot2::ggproto_parent(ggplot2::FacetGrid, self)$draw_panels(...)

    # The finished theme, which is the whole reason this happens here rather
    # than at draw time. See aev_mark_scale().
    theme <- aev_facet_arguments(...)$theme

    # READ BY POSITION, from `FacetGrid$draw_panels`'s own formals, so this is
    # the most fragile lookup in the file. Silently it would leave `scale` at 1
    # and every mark at reference size while the type around them grew.
    aev_require_layout(theme, "the theme ggplot2 hands to a facet")

    # NOT AN ERROR AND NOTHING IS OVERRIDDEN. A user is entitled to restyle a
    # ggplot, and `+` winning is the whole contract; nothing a theme does makes
    # this chart wrong, only off-design. So it says so once and draws anyway.
    aev_warn_if_theme_replaced(theme)

    scale <- aev_mark_scale(theme)

    table <- aev_stamp_scale(table, scale)
    table <- aev_set_panel_heights(table, scale)
    table <- aev_drop_lowest_x_axis(table)
    # After the lower axes are gone, so only the one that is drawn is wrapped.
    table <- aev_wrap_axis_labels(table, self$label_lines)
    # A column variable means groups, and groups mean strips. Without them the
    # brackets and the gap would simply not appear, which is what happened when
    # `residuals = FALSE` and nobody noticed until a review.
    if (length(self$params$cols) > 0L) {
      aev_require_layout(
        which(startsWith(table$layout$name, "strip-b")),
        "the group name strips (`strip-b...`)"
      )
    }

    table <- aev_move_group_strips(table, scale)
    table <- aev_require_panel_titles(table)
    aev_shrink_lower_panel(table)
  }
)

# The rows of the assembled table that hold panels, top to bottom.
aev_panel_rows <- function(table) {

  rows <- sort(unique(table$layout$t[startsWith(table$layout$name, "panel")]))

  # There is always at least one panel. None means the cell naming has changed
  # underneath us, and everything downstream reads the wrong rows.
  aev_require_layout(rows, "any panel (cells named `panel...`)")

  rows
}

# Leaves the x axis between the panels and takes away the one below them, by
# emptying its cell and closing the row up. With a single panel there is
# nothing below to remove and the axis stays where it was.
aev_drop_lowest_x_axis <- function(table) {

  rows <- aev_panel_rows(table)
  if (length(rows) < 2L) return(table)

  below <- which(startsWith(table$layout$name, "axis-b") & table$layout$t > max(rows))

  # `axes = "all_x"` puts one under EVERY panel, so with two panels there is
  # always one below the lower. Without it the chart would print the record
  # names twice and say nothing.
  aev_require_layout(below, "the x axis below the residual panel (`axis-b...`)")

  for (cell in below) {
    table$heights[table$layout$t[cell]] <- grid::unit(0, "cm")
    table$grobs[[cell]] <- ggplot2::zeroGrob()
  }

  table
}

# The row of the x axis that sits between the panels, which is where the group
# names belong too.
aev_interior_axis_row <- function(table) {

  rows <- aev_panel_rows(table)
  if (length(rows) < 2L) return(NA_integer_)

  interior <- table$layout$t[
    startsWith(table$layout$name, "axis-b") &
      table$layout$t > min(rows) &
      table$layout$t < max(rows)
  ]

  if (length(interior) == 0L) NA_integer_ else max(interior)
}

#### Wrapping the record names ####

# WHY THIS IS SPLIT ACROSS BUILD TIME AND DRAW TIME, which is not where anyone
# would choose to put it.
#
# ggplot2 does not wrap axis text at all, and the obvious fix -- a grob that
# breaks its labels to the column width when it is drawn -- founders on the
# layout. Measured on 2026-08-25:
#
# * Inside the x-axis cell, `npc` IS the panel width. The cell spans exactly the
#   panel's columns, so a grob drawn there can break its text to the real width.
# * But `heightDetails`, which is what gtable asks to size the axis row, is
#   called in the enclosing viewport -- the whole device. On an 8 inch chart it
#   sees 8 inches where the panel is 6.67.
# * And a grob that fills itself in at draw time reports a `grobHeight` of zero,
#   so its row would be given no height at all.
#
# So the HEIGHT has to be settled before the width is known, and the BREAKS
# cannot be until after. The split follows: the number of lines is fixed at
# build time, where it becomes real newlines in the labels and ggplot2 measures
# the row correctly; the positions of the breaks are recomputed at draw time
# from the drawn width of the words, into that same number of lines.
#
# THE ONE ASSUMPTION, and it is only ever used to count lines. Measured by
# rendering: the x-axis font is 7.2 pt, a lower-case character averages 0.0476
# inches, and a chart drawn at `ggsave`'s default 7 inches has a panel 5.67
# inches wide. That is 119 characters across the axis. A chart drawn wider
# reserves a line it does not need; one drawn narrower packs the words tighter
# than they want to go, which is what happened to every label before this
# existed.
aev_axis_assumed_chars <- 119

# Past this the axis starts eating the panel it belongs to.
aev_axis_max_lines <- 3L

aev_words <- function(label) {
  words <- strsplit(label, "[[:space:]]+")[[1]]
  words[nzchar(words)]
}

# First-fit packing: take words until the next one would not fit, then break.
# `width_of` is `nchar` at build time and the drawn width at draw time.
aev_pack <- function(words, limit, width_of) {

  lines <- character(0)
  current <- ""

  for (word in words) {
    trial <- if (nzchar(current)) paste(current, word) else word
    if (nzchar(current) && width_of(trial) > limit) {
      lines <- c(lines, current)
      current <- word
    } else {
      current <- trial
    }
  }

  c(lines, current)
}

# The most even break into at most `lines` lines, found by bisecting on the
# width allowed. Packing is monotone in that width -- a wider limit never needs
# more lines -- so the narrowest limit that still fits is the balanced one, and
# bisection finds it.
# As many lines as the name needs at this width, up to the reserved ceiling,
# and the most even break into the ceiling if it needs more than that.
aev_wrap_to_width <- function(label, width, lines, width_of) {

  words <- aev_words(label)
  if (lines <= 1L || length(words) <= 1L) return(label)

  packed <- aev_pack(words, width, width_of)

  if (length(packed) <= lines) return(paste(packed, collapse = "\n"))

  aev_wrap_to_lines(label, lines, width_of)
}

aev_wrap_to_lines <- function(label, lines, width_of) {

  words <- aev_words(label)
  if (lines <= 1L || length(words) <= 1L) return(label)

  low <- max(vapply(words, width_of, numeric(1)))
  high <- width_of(paste(words, collapse = " "))

  for (step in seq_len(40L)) {
    middle <- (low + high) / 2
    if (length(aev_pack(words, middle, width_of)) <= lines) high <- middle else low <- middle
  }

  paste(aev_pack(words, high, width_of), collapse = "\n")
}

# How many lines to reserve, at the assumed width. One number for the whole
# chart: every record is the same width whichever group it is in, and a band
# that changed depth from column to column would look like a mistake.
aev_label_lines <- function(labels) {

  if (length(labels) == 0L) return(1L)

  limit <- aev_axis_assumed_chars / length(labels)
  wanted <- vapply(
    labels,
    function(label) length(aev_pack(aev_words(label), limit, nchar)),
    integer(1)
  )

  min(aev_axis_max_lines, max(1L, max(wanted)))
}

aev_wrap_labels <- function(labels, lines) {
  vapply(labels, aev_wrap_to_lines, character(1),
         lines = lines, width_of = nchar, USE.NAMES = FALSE)
}

# Put a grob back where it came from, found by name. `editGrob` would do this
# through a `gPath`, but the names here are ggplot2's own and auto-generated,
# and matching them by path depends on how strictly grid reads a one-element
# path. Walking is shorter than finding out.
aev_relabel <- function(grob, name, labels) {

  if (identical(grob$name, name)) {
    grob$label <- labels
    return(grob)
  }

  for (slot in c("children", "grobs")) {
    if (!is.null(grob[[slot]])) {
      for (index in seq_along(grob[[slot]])) {
        grob[[slot]][[index]] <- aev_relabel(grob[[slot]][[index]], name, labels)
      }
    }
  }

  grob
}

aev_axis_labels <- function(axis, lines) {
  grid::gTree(children = grid::gList(axis), lines = lines, cl = "aev_axis_labels")
}

#' @exportS3Method grid::makeContent
makeContent.aev_axis_labels <- function(x) {

  axis <- x$children[[1]]
  text <- aev_text_within(axis)

  if (is.null(text) || x$lines <= 1L) {
    return(grid::setChildren(x, grid::gList(axis)))
  }

  labels <- as.character(text$label)

  # THE COLUMN'S REAL WIDTH IS WHAT NAMES ARE PACKED AGAINST. Reserving two
  # lines for the chart says the band is deep enough for two; it does not say
  # every name wants two. Breaking each name into the reserved number regardless
  # split "male lives" across two lines on a chart with room for twenty
  # characters, because one long name elsewhere had forced the reservation.
  # Found by review on 2026-08-26.
  #
  # The reserved count is still the ceiling: a name that will not fit the column
  # in that many lines is broken as evenly as it can be instead, which is what
  # the band was made deep enough for.
  #
  # THE GROB'S OWN `gp`, not the current one: the theme's size and family live
  # on the text and nowhere else, and measuring without them measures a
  # different font.
  width_of <- function(string) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(string, gp = text$gp)), "inches",
      valueOnly = TRUE
    )
  }

  # One record's width. The axis spans its panel's column, and the labels in it
  # are one per record, evenly spaced.
  column <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  record <- column / max(1L, length(labels))

  # From the unwrapped name: the build-time breaks are in the label already and
  # would otherwise be taken as words.
  wrapped <- vapply(
    gsub("\n", " ", labels, fixed = TRUE), aev_wrap_to_width, character(1),
    width = record, lines = x$lines, width_of = width_of, USE.NAMES = FALSE
  )

  grid::setChildren(x, grid::gList(aev_relabel(axis, text$name, wrapped)))
}

# The x axes that are still drawn, which after `aev_drop_lowest_x_axis` is the
# row between the panels.
aev_wrap_axis_labels <- function(table, lines) {

  if (is.na(lines) || lines <= 1L) return(table)

  for (cell in which(startsWith(table$layout$name, "axis-b"))) {
    if (!inherits(table$grobs[[cell]], "zeroGrob")) {
      table$grobs[[cell]] <- aev_axis_labels(table$grobs[[cell]], lines)
    }
  }

  table
}

# THE BRACKET ROUND A GROUP NAME, from `XAxisGroups.Draw` in the C#, which is
# the only place the rule is written down. Reading it off `Test AE chart.svg`
# alone would have given the shape but not the arithmetic.
#
# For each group: a vertical rule at each end of its span, running the full
# depth of the name band, and a horizontal rule across the middle BROKEN either
# side of the name. The break is the name's own width plus a margin at each end,
# and the two remaining segments share what is left equally. Where the name is
# too wide to leave any, the horizontal is dropped and the two verticals carry
# the grouping on their own.
#
# Everything here is a length in drawn text, so it can only be settled once the
# device is known -- which is what `makeContent` is for. The strip that ggplot2
# built is kept as a child rather than redrawn, so the name stays in whatever
# the theme says.
#
# The verticals sit on the edges of the cell, and the cell is the group's own
# panel column: `panel.spacing` is what separates one group from the next, so
# ggplot2 has already narrowed the span the way `widthAdj` does in the C#.
# THE STRIP IS A CHILD FROM THE START, not a field resolved into one. Anything
# that walks a grob for its text -- `labels_in` in the tests, and whatever a
# user reaches in with -- should find the group name in an unresolved bracket
# just as it did before there was a bracket at all.
aev_group_bracket <- function(strip, scale = 1) {
  grid::gTree(children = grid::gList(strip), scale = scale, cl = "aev_group_bracket")
}

# The first piece of drawn text anywhere in a grob. A strip is a `gtable` round
# a `titleGrob` round the text, and the depth is not worth knowing.
aev_text_within <- function(grob) {

  if (inherits(grob, "text")) return(grob)

  for (child in c(grob$children, grob$grobs)) {
    found <- aev_text_within(child)
    if (!is.null(found)) return(found)
  }

  NULL
}

#' @exportS3Method grid::makeContent
makeContent.aev_group_bracket <- function(x) {

  # The bracket is type-sized like everything else, and the strip inside it
  # carries the theme's own font, so a bracket drawn at one weight round a name
  # drawn at another would come apart as the chart grew.
  geometry <- aev_scaled(
    aev_geometry, x$scale %||% 1, c("grid_width_px", "group_gap_px")
  )

  rule <- grid::gpar(
    col = aev_palette$grid_line,
    lwd = geometry$grid_width_px,
    lineend = "butt"
  )

  strip <- x$children[[1]]
  text <- aev_text_within(strip)

  # `grobWidth` rather than `stringWidth`: the strip carries the theme's font
  # in its own `gp`, and only the grob knows it.
  label <- if (is.null(text)) {
    0
  } else {
    grid::convertWidth(grid::grobWidth(text), "npc", valueOnly = TRUE)
  }

  margin <- grid::convertWidth(
    inches(geometry$group_gap_px), "npc", valueOnly = TRUE
  )

  spare <- 1 - (label + 2 * margin)

  ends <- grid::segmentsGrob(
    x0 = grid::unit(c(0, 1), "npc"),
    x1 = grid::unit(c(0, 1), "npc"),
    y0 = grid::unit(c(0, 0), "npc"),
    y1 = grid::unit(c(1, 1), "npc"),
    gp = rule
  )

  across <- if (spare > 0) {
    grid::segmentsGrob(
      x0 = grid::unit(c(0, 1 - spare / 2), "npc"),
      x1 = grid::unit(c(spare / 2, 1), "npc"),
      y0 = grid::unit(c(0.5, 0.5), "npc"),
      y1 = grid::unit(c(0.5, 0.5), "npc"),
      gp = rule
    )
  } else {
    grid::nullGrob()
  }

  grid::setChildren(x, grid::gList(ends, across, strip))
}

# Group names go directly under the record names, which is slide 2's stack:
# plot area, x scale, group names, residual panel.
#
# THE FACET PUTS THEM AT THE VERY BOTTOM and there is no argument that says
# otherwise -- column strips belong to the whole grid, not to one row of it. So
# the row is made where it is wanted and the strips are moved into it, leaving
# the original row empty and closed up.
# WHERE THE GROUP BAND GOES. Between the panels when there are two, and under
# the only x axis when `residuals = FALSE` leaves one -- which used to return NA
# and abandon the strips where the facet put them, so a grouped chart drawn
# without the residual panel had no brackets and no gap at all. Found by review
# on 2026-08-26.
aev_group_destination <- function(table) {

  interior <- aev_interior_axis_row(table)
  if (!is.na(interior)) return(interior)

  axes <- table$layout$t[startsWith(table$layout$name, "axis-b")]
  if (length(axes) == 0L) return(NA_integer_)

  max(axes)
}

aev_move_group_strips <- function(table, scale = 1) {

  strips <- which(startsWith(table$layout$name, "strip-b"))
  if (length(strips) == 0L) return(table)

  # There are strips to place, so there must be a row to place them in.
  destination <- aev_group_destination(table)
  aev_require_layout(destination, "a row for the group names (`axis-b...`)")

  # Captured before anything moves, along with the columns they span.
  moving <- table$grobs[strips]
  left <- table$layout$l[strips]
  right <- table$layout$r[strips]

  vacated <- unique(table$layout$t[strips])
  height <- table$heights[vacated[1]]

  for (cell in strips) table$grobs[[cell]] <- ggplot2::zeroGrob()
  table$heights[vacated] <- grid::unit(0, "cm")

  # Two rows: the gap, then the band. See `group_band_gap_px`.
  table <- gtable::gtable_add_rows(
    table, inches(aev_geometry$group_band_gap_px * scale), pos = destination
  )
  table <- gtable::gtable_add_rows(table, height, pos = destination + 1L)

  for (k in seq_along(moving)) {
    table <- gtable::gtable_add_grob(
      table, aev_group_bracket(moving[[k]], scale),
      t = destination + 2L, l = left[k], r = right[k],
      name = paste0("aev-group-", k), clip = "off"
    )
  }

  table
}

# NOTHING IS DONE TO THE A/E PANEL'S TITLE ANY MORE. The theme carries its full
# size -- `YAxisTitleFontSize`, 16 px -- where it used to be sized like an axis
# number and scaled up here afterwards. What is still worth doing is insisting
# that it EXISTS: a chart quietly missing a panel name is the sort of wrong that
# looks like a design choice.
aev_require_panel_titles <- function(table) {

  rows <- aev_panel_rows(table)

  titles <- which(
    startsWith(table$layout$name, "strip-l") & table$layout$t == min(rows)
  )
  aev_require_layout(titles, "the A/E panel's title (`strip-l...`)")

  table
}

# THE RESIDUAL PANEL'S DEPTH FOLLOWS THE TYPE, like everything else on the
# chart. Left fixed, the strip stays 84 px however large the type is, so a chart
# asking for half again the type gets the same seven rows of residual in the
# same depth with everything in them half again as big -- which reads as the
# residual axis tightening up while the A/E panel above it opens out. Tim spotted
# exactly that on 2026-08-25.
#
# WHY THE HEIGHTS ARE SET HERE AND NOT BY `theme(panel.heights = )`, which is
# where they used to be set and is the obvious place: ggplot2 applies that theme
# setting AFTER the facet has drawn, so anything the facet does to those rows is
# overwritten. Measured, not guessed -- the depth stayed at 84 px. Setting them
# here instead leaves nothing to overwrite them, and only the facet has the
# scale anyway, since the theme's own `base_size` misses a later `+ theme(...)`.
aev_set_panel_heights <- function(table, scale) {

  rows <- aev_panel_rows(table)
  if (length(rows) == 0L) return(table)

  heights <- aev_panel_heights(length(rows) >= 2L, scale)

  for (index in seq_along(rows)) {
    table$heights[rows[index]] <- heights[index]
  }

  table
}

# The residual panel's numbers and title, smaller than the A/E panel's.
#
# Found by row rather than by name: whatever sits to the left of the bottom
# panel is that panel's scale and title, however the cells are numbered.
aev_shrink_lower_panel <- function(table) {

  rows <- aev_panel_rows(table)
  if (length(rows) < 2L) return(table)

  lowest <- max(rows)
  in_lowest <- function(prefix) {
    which(startsWith(table$layout$name, prefix) & table$layout$t == lowest)
  }

  # BOTH THE NUMBERS AND THE TITLE, each by its own ratio. The theme carries
  # the A/E panel's sizes; this panel differs from it in two ways and they are
  # not the same factor -- `scale_ratio` makes the numbers BIGGER (Tim, on the
  # grounds that this axis is always the same seven values) while
  # `title_ratio` makes the title smaller, which is `DevYAxisTitleRatio`.
  title <- in_lowest("strip-l")
  aev_require_layout(title, "the residual panel's title (`strip-l...`)")

  table$grobs[title] <- lapply(
    table$grobs[title], aev_scale_text, aev_residual$title_ratio
  )

  scale <- in_lowest("axis-l")

  # The residual panel always has numbers down its side, and they are always
  # smaller than the A/E panel's. Missing them leaves the two panels' scales
  # the same size, which reads as a chart that was simply designed that way.
  aev_require_layout(scale, "the residual panel's axis numbers (`axis-l...`)")

  table$grobs[scale] <- lapply(table$grobs[scale], aev_scale_text, aev_residual$scale_ratio)

  table
}

# `cex` rather than `fontsize`, because it multiplies whatever size each piece
# of text already has instead of overriding it, and it is inherited down a
# tree. Merged into any existing settings rather than replacing them.
aev_scale_text <- function(grob, ratio) {
  gp <- grob$gp %||% grid::gpar()
  gp$cex <- ratio
  grob$gp <- gp
  grob
}

#### The legend ####

# WHAT THE FOUR ENTRIES ARE FOR. Three of them take the marker apart -- the
# point, the outer bar and the inner one -- because the marker is three
# statements in one shape and nothing else on the chart says which is which.
# The fourth names the panel below.
#
# `Legend0..3` in `AEChart.cs`, in that order, and the same order in the
# rendered samples. The gloss on the residual is NOT here: Tim put it in the
# panel title, where it sits against the thing it explains.
aev_legend_entries <- function(residuals) {

  # A LIST, NOT A CHARACTER VECTOR, because one entry has to be plotmath and the
  # others must not be. Sigma is U+03C3, which is not in latin1 and so cannot be
  # drawn from a string on the device `R CMD check` uses -- see
  # `aev_signed_math()` for the mechanism. `expression(sigma)` takes it from the
  # Symbol font instead, which is how base R has always drawn Greek.
  entries <- list(
    "log A/E",
    "95% confidence",
    expression(paste(pm, sigma, " (~68%) confidence")),
    "Deviance residual"
  )

  if (residuals) entries else entries[seq_len(3)]
}

# ONE LAYER, FOUR GLYPHS. ggplot2 draws every layer's key into every entry of a
# legend, so four layers would give four entries each showing all four glyphs
# on top of each other. A single layer whose `draw_key` looks at which entry it
# has been handed is the way round it, and the mapped `shape` is what says
# which -- it is never drawn, only read.
#
# The layer itself draws nothing at all: its x and y are missing. It exists so
# that a scale exists, so that a legend exists.
aev_legend_layer <- function(residuals, panel, run) {

  entries <- aev_legend_entries(residuals)

  data <- data.frame(
    x = rep(NA_real_, length(entries)),
    y = NA_real_,
    key = factor(entries, levels = entries),
    panel = panel
  )
  data$run <- run

  ggplot2::geom_point(
    data = data,
    mapping = ggplot2::aes(x = .data$x, y = .data$y, shape = .data$key),
    na.rm = TRUE,
    key_glyph = aev_draw_key
  )
}

# THE KEY IS A SAMPLE OF THE CHART, so it has to be drawn at the chart's size.
#
# The legend is the one part of the plot the facet never sees -- it is added
# after `draw_panels` has returned -- so the scale cannot be stamped on it the
# way it is stamped on the marks. What the key IS told is its own box, in
# millimetres, in `size`. Measured on ggplot2 4.0.3: 6.35 x 4.23 mm by default,
# and it follows `legend.key.width` and `legend.key.height` from any source.
#
# So the glyph is sized against its box, and the box is sized against the type
# in `theme_aev()`. The fontsize in force at this point is the device's 12
# whatever the theme says -- checked, not assumed -- so it is no use here.
aev_key_scale <- function(size) {

  reference <- aev_legend$key_width_px / 96 * 25.4

  if (length(size) == 0L) return(1)

  width <- as.numeric(size[1])
  if (!is.finite(width) || width <= 0) return(1)

  width / reference
}

aev_draw_key <- function(data, params, size) {

  scale <- aev_key_scale(size)

  switch(
    data$shape[1],
    aev_key_marker(scale),
    aev_key_confidence(scale),
    aev_key_sigma(scale),
    aev_key_residual(scale)
  )
}

# The glyphs, from the four `<line>`s and two `<circle>`s that open
# `Test AE chart.svg`. Each is the piece of the marker it names, at the size it
# is drawn on the chart -- so the key is a sample rather than a diagram.
#
# Across the key the positions are proportions and the sizes are pixels, which
# is the same division the marker itself makes.
aev_key_middle <- grid::unit(0.5, "npc")

aev_key_marker <- function(scale = 1) {

  geometry <- aev_scaled(
    aev_geometry, scale, c("marker_radius_px", "thin_width_px")
  )
  radius <- inches(geometry$marker_radius_px)

  grid::grobTree(
    grid::segmentsGrob(
      x0 = aev_key_middle, x1 = aev_key_middle,
      y0 = grid::unit(c(0, 1), "npc"),
      y1 = grid::unit.c(aev_key_middle - radius, aev_key_middle + radius),
      gp = grid::gpar(
        col = aev_palette$marker_weak,
        lwd = geometry$thin_width_px,
        lineend = "butt"
      )
    ),
    grid::circleGrob(
      x = aev_key_middle, y = aev_key_middle, r = radius,
      gp = grid::gpar(
        col = aev_palette$marker_weak,
        fill = NA,
        lwd = geometry$thin_width_px
      )
    )
  )
}

aev_key_confidence <- function(scale = 1) {

  geometry <- aev_scaled(
    aev_geometry, scale, c("cap_width_px", "cap_height_px", "thick_width_px")
  )
  half_cap <- inches(geometry$cap_width_px / 2)
  cap <- grid::unit(0.7, "npc")

  grid::grobTree(
    grid::segmentsGrob(
      x0 = aev_key_middle, x1 = aev_key_middle,
      y0 = cap, y1 = grid::unit(0.15, "npc"),
      gp = grid::gpar(
        col = aev_palette$marker_medium,
        lwd = geometry$thick_width_px,
        lineend = "butt"
      )
    ),
    grid::segmentsGrob(
      x0 = aev_key_middle - half_cap, x1 = aev_key_middle + half_cap,
      y0 = cap, y1 = cap,
      gp = grid::gpar(
        col = aev_palette$marker_cap,
        lwd = geometry$cap_height_px,
        lineend = "butt"
      )
    )
  )
}

aev_key_sigma <- function(scale = 1) {

  geometry <- aev_scaled(aev_geometry, scale, c("thin_width_px", "thick_width_px"))


  grid::grobTree(
    grid::segmentsGrob(
      x0 = aev_key_middle, x1 = aev_key_middle,
      y0 = grid::unit(0.95, "npc"), y1 = aev_key_middle,
      gp = grid::gpar(
        col = aev_palette$marker_medium,
        lwd = geometry$thick_width_px,
        lineend = "butt"
      )
    ),
    grid::segmentsGrob(
      x0 = aev_key_middle, x1 = aev_key_middle,
      y0 = aev_key_middle, y1 = grid::unit(0.05, "npc"),
      gp = grid::gpar(
        col = aev_palette$marker_weak,
        lwd = geometry$thin_width_px,
        lineend = "butt"
      )
    )
  )
}

aev_key_residual <- function(scale = 1) {

  residual <- aev_scaled(
    aev_residual, scale,
    c("marker_radius_px", "marker_width_px", "drop_width_px")
  )

  radius <- inches(residual$marker_radius_px)
  centre <- grid::unit(0.72, "npc")

  grid::grobTree(
    grid::segmentsGrob(
      x0 = aev_key_middle, x1 = aev_key_middle,
      y0 = centre - radius, y1 = grid::unit(0.1, "npc"),
      gp = grid::gpar(
        col = aev_palette$marker_medium,
        lwd = residual$drop_width_px[1],
        lineend = "butt"
      )
    ),
    grid::circleGrob(
      x = aev_key_middle, y = centre, r = radius,
      gp = grid::gpar(
        col = aev_palette$marker_medium,
        fill = NA,
        lwd = residual$marker_width_px
      )
    )
  )
}

# Slide 2's layout rules, as far as a theme can carry them:
#
# * Text is never rotated, so the panel titles lie flat and are left-aligned.
# * Plot lines, ticks included, never appear outside the plot area, so there
#   are no tick marks.
# * X-axis grid lines are shown only for non-categorical charts, and this one
#   is categorical.
#
# THE Y-AXIS TITLES ARE FACET STRIPS. Two panels need two titles and a plot has
# only one `axis.title.y`, so the panel name is the title and `switch = "y"`
# puts it where slide 2's layout stack puts it. `strip.placement = "outside"`
# keeps it clear of the scale.
#
# THE GRID LINES ARE NOT IN THE THEME. They are drawn as a layer restricted to
# the A/E panel, because a theme grid line would appear in the residual panel
# too, where `AEChart.cs` deliberately has none -- the tint band edges do that
# work, and they fall in the same places.
# The name of the element the chart signs its theme with. Registered in
# `.onLoad` -- see `R/zzz.R` for why it is an element rather than an attribute.
aev_signature <- "logmu.signature"

# HAS THE CHART'S OWN THEME BEEN REPLACED? Only a complete theme clears the
# signature, so `+ theme(...)` never trips this and `+ theme_minimal()` always
# does.
#
# THE TEST IS "NOT TRUE", NOT "IS NULL", and the difference is the whole reason
# this took two goes. A raw replaced theme reads NULL, but the theme that
# reaches the facet has been through `plot_theme()`, which completes it against
# the registered defaults -- so the signature comes back as the registered
# FALSE, not as NULL, and a test for NULL never fired.
#
# A missing registration is treated as "cannot tell" rather than as a
# replacement: warning because the package failed to load properly would be
# worse than saying nothing.
aev_theme_replaced <- function(theme) {

  # NO THEME IS NOT A REPLACED THEME. `aev_mark_scale()` tolerates NULL the same
  # way and for the same reason: if the theme cannot be read at all, that is a
  # fault to stay quiet about rather than to blame on the user's styling.
  if (is.null(theme)) return(FALSE)

  !isTRUE(tryCatch(
    ggplot2::calc_element(aev_signature, theme),
    error = function(condition) TRUE
  ))
}

aev_warn_if_theme_replaced <- function(theme) {

  if (!aev_theme_replaced(theme)) return(invisible(FALSE))

  warning(
    "This chart's theme has been replaced, so some of its design is lost: ",
    "grid lines return in the residual panel, where the tint bands already ",
    "mark the same places, and the panel titles are rotated upright instead ",
    "of lying flat beside the panel. Add `theme_aev()` after the other theme ",
    "to restore it.",
    call. = FALSE
  )

  invisible(TRUE)
}

#' The chart theme, and how to resize a chart
#'
#' @description
#' The theme [autoplot.aev()] draws with. Adding it back at a different
#' `base_size` is how a **logmu** chart is resized.
#'
#' @details
#' Every length on the chart is a multiple of the base text size -- the marker
#' and its interval, the off-scale chevrons, the depth of the residual strip,
#' the group brackets and the legend key. So a chart drawn for a larger page
#' wants larger type, and gets everything else with it:
#'
#' ```
#' autoplot(aev) + theme_aev(base_size = 13)
#' ```
#'
#' A larger device on its own does not do this. It gives the same marker in more
#' room, which reads thinner and smaller the further it is stretched.
#'
#' `ggplot2::theme(text = ggplot2::element_text(size = 13))` resizes most of the
#' chart in the same way and is the ordinary **ggplot2** way of asking. It
#' cannot reach the legend key, which is sized by `legend.key.width` and
#' `legend.key.height` and so only by this function. Where the legend matters,
#' use `theme_aev()`.
#'
#' The result is an ordinary theme and can be extended with `+` like any other.
#'
#' @section Adding a standard theme:
#' Replacing this one wholesale -- with `ggplot2::theme_minimal()`, say -- still
#' draws a chart, and rather more of one than you might expect. The two panels,
#' the group brackets, the off-scale chevrons, the depth of the residual strip
#' and the sizing all survive, because the facet builds them and not the theme.
#'
#' What is lost is the styling the design asks for:
#'
#' * Grid lines come back in the residual panel, where the prototype
#'   deliberately has none -- the edges of the tint bands do that work, and they
#'   fall in the same places.
#' * `ggplot2::theme_bw()` and `ggplot2::theme_grey()` restore the tick marks
#'   too.
#' * The panel titles are facet strips, so they revert to being rotated upright
#'   and centred inside the panel instead of lying flat beside it.
#' * The type goes from 9 points to 11, which takes every mark on the chart with
#'   it -- markers about a fifth heavier, and a residual strip a fifth deeper.
#'
#' None of that is an error, and nothing is overridden: a chart is entitled to
#' be restyled, and `+` winning is the whole of **ggplot2**'s contract. But the
#' chart does say so. It signs its own theme, notices when the signature has
#' gone, and warns once as it draws:
#'
#' ```
#' This chart's theme has been replaced, so some of its design is lost ...
#' ```
#'
#' Only a COMPLETE theme clears the signature. `ggplot2::theme(...)` merges, so
#' ordinary customising never trips it. Adding `theme_aev()` after the standard
#' theme puts the styling back and silences the warning.
#'
#' @param base_size Base text size in points. Every other size on the chart is a
#'   multiple of it. The default of 9 is the size the design was drawn at, and
#'   at that size every length is the pixel the prototype specifies.
#' @returns A **ggplot2** theme.
#' @seealso [aev_plot] for the chart itself.
#' @examples
#' aev <- create_aev(
#'   A = c(1100, 40, 260),
#'   E = c(1000, 50, 250),
#'   V = c(2500, 125, 620)
#' )
#' names(aev) <- c("[65, 75)", "[75, 85)", "[85, 95)")
#'
#' # The same chart, drawn for a larger page.
#' autoplot(aev) + theme_aev(base_size = 13)
#' @export
theme_aev <- function(base_size = 9) {

  type_scale <- base_size / aev_reference_text_size

  # EVERY SIZE ON THIS CHART IS THE DESIGN'S OWN, AS A RATIO OF THE BASE.
  #
  # `ChartStyle` states its sizes in CSS pixels against a 12 px base, and this
  # chart's base text is 9 pt, which IS that 12 px. So the ratio is simply the
  # design's pixels over twelve, and every one of them follows `base_size`.
  #
  # `rel()` COMPOSES DOWN THE INHERITANCE CHAIN rather than replacing what it
  # inherits. `theme_minimal()` sets `strip.text`, `axis.text` and `legend.text`
  # to rel(0.8), so without neutralising those three first, every ratio below
  # would quietly come out at four fifths of what it says -- which is exactly
  # how the group names ended up at 8.28 pt while claiming 1.15 of the base.
  design <- function(pixels) ggplot2::rel(pixels / 12)

  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      strip.text  = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text   = ggplot2::element_text(size = ggplot2::rel(1))
    ) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = aev_palette$plot_background, colour = NA),
      plot.title = ggplot2::element_text(
        colour = aev_palette$title_text,
        hjust = 0,
        vjust = 0,
        size = design(aev_headings$title_px),
        margin = ggplot2::margin(
          b = base_size * aev_headings$margin_ratio, unit = "pt"
        )
      ),
      plot.title.position = "plot",
      # A caption is where the overdispersion goes, since an `aev` does not
      # carry it and the chart cannot know it. Left-aligned with the title
      # rather than in ggplot2's bottom-right, because slide 2 left-aligns
      # everything of this kind.
      plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(
        colour = aev_palette$axis_scale_text,
        hjust = 0,
        size = design(aev_notes$text_px),
        lineheight = aev_notes$line_height,
        # A margin rather than a blank first note: the gap belongs between the
        # chart and the block, not inside the text.
        margin = ggplot2::margin(t = base_size * aev_notes$gap_ratio, unit = "pt")
      ),
      # A SUBTITLE IS STYLED EVEN THOUGH NOTHING SETS ONE BY DEFAULT. It works
      # through `subtitle`, and without this it came out in plain black while
      # every other piece of text on the chart is off the palette.
      plot.subtitle = ggplot2::element_text(
        colour = aev_palette$axis_title_text,
        hjust = 0,
        size = design(aev_headings$subtitle_px),
        margin = ggplot2::margin(
          b = base_size * aev_headings$margin_ratio, unit = "pt"
        )
      ),
      panel.background  = ggplot2::element_rect(fill = aev_palette$plot_background, colour = NA),
      # NO FRAME ROUND THE PANELS (Tim, 2026-08-21). The outermost grid lines
      # already close the top and bottom of the A/E panel, and with groups the
      # borders drew a box round every column, which read as a grid of boxes
      # rather than as one chart cut into groups.
      panel.border      = ggplot2::element_blank(),
      panel.grid        = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_blank(),
      # `XAxisScaleFontSize` and `YAxisScaleFontSize`, both 9 px.
      axis.text.x = ggplot2::element_text(
        colour = aev_palette$axis_scale_text, vjust = 1, size = design(9)
      ),
      axis.text.y = ggplot2::element_text(
        colour = aev_palette$axis_scale_text, hjust = 1, size = design(9)
      ),
      axis.title.y      = ggplot2::element_blank(),
      strip.placement   = "outside",
      strip.background  = ggplot2::element_blank(),
      # `YAxisTitleFontSize`, 16 px. These are the panel titles, and they are
      # axis titles rather than axis numbers -- the design sizes them as such,
      # and sizing them like numbers was what made them read as an afterthought.
      strip.text.y.left = ggplot2::element_text(
        colour = aev_palette$axis_title_text,
        angle = 0,
        hjust = 0,
        size = design(16)
      ),
      # Slide 5 says group names take "the same font size and colour as the
      # y-axis title", which is the strip above. TIM OVERRODE THE COLOUR ON
      # 2026-08-25: a group name labels the records under it, so it belongs to
      # the same voice as their names rather than to the axis titles.
      # `XAxisGroupNameFontSize`, 12 px. That is a third larger than the record
      # names at 9, which is the distinction Tim asked for on 2026-08-25 and
      # more of it than the 1.15 he was shown.
      strip.text.x.bottom = ggplot2::element_text(
        colour = aev_palette$axis_scale_text,
        size = design(12),
        angle = 0
      ),
      panel.spacing.x = inches(aev_geometry$group_gap_px),
      # Slide 2 puts the legend above the plot area and below the title, and
      # left-aligns its items with the other two. Slide 5 adds that the "key is
      # small and faded", which is the size and the colour here.
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "left",
      # "panel", not "plot": the legend lines up with the y axis rather than
      # with the far left edge, which is where Tim wants its left edge.
      legend.location = "panel",
      legend.title = ggplot2::element_blank(),
      # A box round the legend, in the weight the panels used to have. The
      # slides show none and the C# has `LegendBorderWidth`; Tim chose the box
      # on 2026-08-21, which is also what sample B does.
      legend.background = ggplot2::element_rect(
        fill = NA,
        colour = aev_palette$grid_line,
        linewidth = line_width(aev_geometry$border_width_px)
      ),
      legend.key = ggplot2::element_blank(),
      # THE BOX THE KEY IS DRAWN IN, sized against the type like everything
      # else. `aev_draw_key()` reads its own box back and draws the glyph to
      # match, so this one setting carries the whole legend.
      legend.key.width = inches(aev_legend$key_width_px * type_scale),
      legend.key.height = inches(aev_legend$key_height_px * type_scale),
      # `LegendFontSize`, 10 px.
      legend.text = ggplot2::element_text(
        colour = aev_palette$marker_medium,
        size = design(10)
      ),
      # THE SIGNATURE. See `R/zzz.R`. It rides on this theme and is lost only
      # when the whole theme is replaced, which is exactly what the chart wants
      # to be told about.
      logmu.signature = TRUE
    )
}

# The rows of the residual panel that carry a tint, as `Background1To2` and
# `Background2To3`, mirrored above and below the line.
aev_residual_bands <- function(panel) {
  alpha <- aev_residual$band_alpha
  data.frame(
    panel = panel,
    lower = c(1, 2, -2, -3),
    upper = c(2, 3, -1, -2),
    alpha = c(alpha[1], alpha[2], alpha[1], alpha[2])
  )
}

# ARGUMENTS THAT ARE NOT ARGUMENTS.
#
# `...` is there so `plot()` can forward to `autoplot()` and for no other
# reason. Silently dropping whatever else lands in it means
# `plot(aev, main = "Male retirees")` -- which is exactly what somebody arriving
# from base R graphics will write -- draws a chart with no title and says
# nothing about why. Same for a misspelling: `notes` for `note`, `residual` for
# `residuals`.
#
# So anything unrecognised is an error, and the message says what the chart
# does take. The base R names are worth translating by hand, because that is
# the mistake most likely to be made and the least likely to be self-evident.
aev_chart_arguments <- c(
  "title", "subtitle", "notes", "residuals", "log_range", "log_step"
)

# ONLY NAMES THAT MATCH NOTHING GET HERE. R's partial matching runs first, so
# `sub` reaches `subtitle` and `note` reaches `notes` without ever being seen as
# unused -- which is ordinary R behaviour and, in `sub`'s case, a happy accident
# for anyone arriving from base graphics. Listing them here would be dead
# entries claiming to help.
aev_argument_advice <- c(
  main    = "title",
  caption = "notes",
  xlab    = "names(aev)",
  ylab    = "nothing -- the panels name themselves",
  cex     = "theme_aev(base_size = )",
  col     = "nothing -- the palette is fixed"
)

aev_check_dots <- function(...) {

  extra <- list(...)
  if (length(extra) == 0L) return(invisible(NULL))

  given <- names(extra)
  if (is.null(given)) given <- rep("", length(extra))

  named <- given[nzchar(given)]
  unnamed <- sum(!nzchar(given))

  advice <- character(0)

  for (name in named) {
    # `[[` on a name the vector does not hold is an error, not a NULL, so the
    # membership test has to come first.
    advice <- c(advice, if (name %in% names(aev_argument_advice)) {
      paste0("`", name, "` (use ", aev_argument_advice[[name]], ")")
    } else {
      paste0("`", name, "`")
    })
  }

  if (unnamed > 0L) {
    advice <- c(advice, paste(unnamed, "unnamed argument(s)"))
  }

  stop(
    "Unused argument(s) to an `aev` chart: ", paste(advice, collapse = ", "), ".\n",
    "It takes ", paste0("`", aev_chart_arguments, "`", collapse = ", "),
    ", and anything else through `+` as an ordinary ggplot2 chart.",
    call. = FALSE
  )
}

#' Chart an `aev`
#'
#' @description
#' Draws A/E for each record of an `aev` on a log scale, with its confidence
#' interval, as a **ggplot2** object that can be extended with `+` like any
#' other.
#'
#' `plot()` is the same chart, drawn immediately.
#'
#' @section The marker:
#' Each record is an open circle at A/E, with the interval drawn around it in
#' three weights: a thin line out to one standard deviation, a thicker line from
#' there to 1.96 standard deviations, and a bar across the end.
#'
#' The interval is \eqn{\sqrt{V} / E}, symmetric in \eqn{\log(A/E)} and so
#' asymmetric in A/E itself. It says how precisely the record was measured. It
#' is not the significance test: read that off the deviance residual, which is
#' printed by [aev_properties] and drawn in the panel below. Reading whether an
#' interval crosses the 100% line is quick and nearly right; the residual is
#' exact.
#'
#' @section The residual panel:
#' Unless `residuals` is `FALSE`, the deviance residual for each record is drawn
#' below the chart against the same x axis. This is the test the marker only
#' approximates. Its axis is fixed at \eqn{\pm 3.5} whatever the data does, so
#' one chart can be read against another, and the tint bands mark 1 to 2 and 2
#' to 3.
#'
#' The drop line thickens and darkens as the residual grows, at breakpoints of
#' 0.5, 1.5 and 2.5, so that a bad cell is visible before any number is read.
#' Above 2.5 nothing further happens: a residual of 6 is drawn exactly like a
#' residual of 3. That is deliberate. Line weight and lightness carry only a
#' few steps a reader can tell apart, and they are better spent on the range
#' where the reading turns from unremarkable to worth investigating than on a
#' tail where the response is the same either way, where the normal
#' approximation is at its worst, and where every value scales as
#' \eqn{1 / \sqrt{\Omega}} in an overdispersion the caller chose.
#'
#' So do not rank records by how heavy their lines look. The residual is
#' printed by every marker to two places, and that is what separates them.
#'
#' The printed value sits across the zero line from its own marker, so it is
#' never clipped. A residual beyond \eqn{\pm 3.5} has its circle cut off at the
#' edge of the panel, with its drop line running to the boundary, and its
#' number is still there to be read.
#'
#' @section Records with nowhere to go:
#' A record is drawn only if it has a position. Four cases do not:
#'
#' * A, E and V all zero, which is a true statement about nobody.
#' * Any of A, E or V missing.
#' * A of zero with exposure, where \eqn{\log(A/E)} is infinite.
#' * A/E give or take one standard deviation lying wholly off the axis, which
#'   by default runs from 61% to 165%.
#'
#' Each keeps its place on the x axis and its label, so no category disappears.
#'
#' The last two are real answers with nowhere to put them, and they are marked
#' with an orange chevron against the top or bottom of the panel saying which
#' way the record went. A record with no deaths gets one at the bottom, since
#' \eqn{\log(A/E)} of \eqn{-\infty} is below any axis. The first two are not
#' marked: a chevron says which way a record went, and neither "nobody" nor
#' "not known" went anywhere.
#'
#' The test is on the **inner** interval, not the 95% one. A record whose marker
#' and inner interval are both off the axis is replaced by a chevron even where
#' its outer arms would have reached back into the panel, because a bar with no
#' marker and no cap says less than the chevron does. A record whose inner
#' interval does reach the panel is drawn as it stands, with whatever runs past
#' the edge clipped. Widening `log_range` therefore removes chevrons, and
#' narrowing it adds them.
#'
#' @section The names below the chart:
#' Record names sit between the two panels, and group names below them inside a
#' bracket: a vertical rule at each end of the group's span and a horizontal
#' rule across the middle, broken either side of the name. Where the name is too
#' wide to leave any horizontal, the two verticals carry the grouping alone.
#'
#' Record names are wrapped onto more than one line when they will not fit the
#' width of a record. Nothing needs to be passed for this and there is nothing
#' to tune. How many lines to allow is decided from the names and the number of
#' records, against a chart the size `ggsave()` draws by default; where the
#' breaks fall is decided when the chart is drawn, from the width of the words
#' in the font they end up in. A chart drawn much narrower than that will pack
#' its names tighter than they want to go, which is what every chart did before.
#'
#' @section Title, subtitle and notes:
#' Three pieces of text can be set on a chart, and they sit where a reader
#' expects them: the `title` above it, the `subtitle` under the title, and
#' `notes` beneath everything. All three are left-aligned with the plot.
#'
#' `notes` is a character vector and nothing more: each element is printed on
#' its own line beneath the chart, in the order given. Notes are set smaller
#' than anything else on the chart and their lines further apart, so that a list
#' of them reads as separate statements rather than as a paragraph.
#'
#' No heading is added. A note is whatever the caller says it is, so a caller
#' who wants one writes it as the first note.
#'
#' ```
#' autoplot(aev, notes = c(
#'   "Overdispersion 2.0",
#'   "Z = 1.31 against S4PMA",
#'   "Experience to 31 December 2025"
#' ))
#' ```
#'
#' Deliberately, it knows nothing about what a note says. An `aev` is three
#' vectors of numbers and anybody may build one, so the chart is not coupled to
#' how those numbers were produced -- and its notes are not either. Whatever a
#' caller wants recorded against the chart, they write and pass.
#'
#' The notes are the plot's caption, so `+ ggplot2::labs(caption = ...)` added
#' afterwards replaces them, as it would any caption.
#'
#' @section The A/E axis:
#' `log_range` and `log_step` are in \eqn{\log(A/E)}, because that is the space
#' the axis is even in: the default range of \eqn{\pm 0.5} in steps of
#' \eqn{0.1} comes out as 61, 67, 74, 82, 90, 100, 111, 122, 135, 149 and 165
#' per cent.
#'
#' `log_step` must divide both ends of `log_range`. The axis has no expansion,
#' so its ends are its outermost grid lines and a step that fell short would
#' leave the panel open at the top; dividing both ends also guarantees a break
#' on 100%, which is the line every marker is read against. A step of 0.2 does
#' not divide the default range, so twenty per cent steps want a range of
#' \eqn{\pm 0.6}. The error says so when it happens.
#'
#' `log_range` must stay within \eqn{\pm 1}, an A/E of 37% to 272%. The two
#' panels share one scale and are told apart by how far each reaches, so an A/E
#' axis cannot be allowed to grow into the residual panel's territory.
#'
#' These arguments exist because the usual way round is a trap: adding
#' `+ ggplot2::scale_y_continuous(...)` replaces the whole scale, and with it
#' the per-panel limits that keep the residual axis at \eqn{\pm 3.5}. That
#' produces a wrong chart rather than an error.
#'
#' @section Its size:
#' Everything drawn on the chart is sized against its type rather than against
#' the device: the marker and its interval, the chevrons, the residual markers,
#' the depth of the residual strip, the group brackets and the legend key are
#' all multiples of the base text size. So the way to draw this chart larger is
#' to ask for larger type as well as a larger device, which is what
#' [theme_aev()] is for:
#'
#' ```
#' autoplot(aev) + theme_aev(base_size = 13)
#' ```
#'
#' A larger device on its own gives the same marker in more room, which reads
#' thinner and smaller the further it is stretched; raising the type with it
#' keeps the chart in the proportions it was drawn in. There is no separate
#' argument for this, deliberately -- a second knob could be left out of step
#' with the theme, and this one cannot.
#'
#' The A/E panel takes whatever height is left once the title, legend, x axis,
#' group names and residual strip have taken theirs, so it is the part that
#' absorbs the slack. It wants room: the prototype runs at about four times the
#' depth of the residual strip, which a default `ggsave()` comfortably gives.
#'
#' One theme setting is a trap, in the same way a replaced scale is. Adding
#' `+ ggplot2::theme(panel.heights = ...)` overrides both panels' heights after
#' the chart has set them, so the residual strip stops being the fixed depth it
#' is meant to be and stops following the type. Nothing warns: the setting is
#' applied later than anything this chart can see. Resize with
#' [theme_aev()] instead.
#'
#' @param object,x An `aev`.
#' @param title A title for the chart, or `NULL` for none.
#' @param subtitle A subtitle, printed under the title, or `NULL` for none.
#' @param notes Lines of text to print under the chart, as a character vector,
#'   or `NULL` for none. See *Notes*.
#' @param residuals Whether to draw the deviance residual panel below the A/E
#'   chart.
#' @param log_range,log_step The A/E axis, in \eqn{\log(A/E)}: the range it
#'   covers and the spacing of its labelled ticks. See *The A/E axis*.
#' @param ... Passed from `plot()` to `autoplot()`. Anything else is an error
#'   rather than being quietly dropped, so that `plot(aev, main = "...")` says
#'   what it should have been.
#' @returns `autoplot()` returns a `ggplot`. `plot()` draws it and returns `x`
#'   invisibly.
#' @examples
#' aev <- create_aev(
#'   A = c(1100, 40, 260),
#'   E = c(1000, 50, 250),
#'   V = c(2500, 125, 620)
#' )
#' names(aev) <- c("[65, 75)", "[75, 85)", "[85, 95)")
#'
#' autoplot(aev, title = "Male retirees, amounts-weighted")
#'
#' # A wider axis in coarser steps, and the assumption the chart was drawn on.
#' autoplot(aev, log_range = c(-0.6, 0.6), log_step = 0.2) +
#'   ggplot2::labs(caption = "Overdispersion 2.0")
#' @name aev_plot
NULL

#' @rdname aev_plot
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @rdname aev_plot
#' @importFrom ggplot2 .data
#' @export
autoplot.aev <- function(object,
                         title = NULL,
                         subtitle = NULL,
                         notes = NULL,
                         residuals = TRUE,
                         log_range = c(-0.5, 0.5),
                         log_step = 0.1,
                         ...) {

  if (!is_aev(object)) {
    stop("The argument `object` must be an `aev`.", call. = FALSE)
  }

  aev_check_dots(...)

  aev_check_log_scale(log_range, log_step)
  log_breaks <- aev_log_breaks(log_range, log_step)

  panels <- if (residuals) aev_panels else aev_panels[1]
  # `rep` rather than a bare value: assigning a length-one factor to a
  # zero-row data frame is an error, and a zero-record aev is legal.
  as_panel <- function(name, times = 1L) factor(rep(name, times), levels = panels)

  # THE AXIS DECIDES WHAT IS OFF IT. A record outside the chosen range is off
  # scale by definition, so the status has to be judged against `log_range` and
  # not against the default one.
  frame <- aev_plot_frame(object, log_range)
  frame$panel <- as_panel(aev_panels[1], nrow(frame))

  # GROUPS IN RECORD ORDER, not alphabetical: a breakdown of age bands followed
  # by cohorts should read in that order, and `factor()` would sort them.
  # A RUN OF RECORDS, NOT A NAME. The chart is split on where the group name
  # CHANGES, which is what the prototype does -- `XAxisGroups.Draw` walks the
  # records in order and closes a run when the next name differs. Splitting on
  # the name itself instead gathers every record that shares a name into one
  # place however far apart they are, and a name that occurs twice then leaves
  # two stretches of the axis overlapping: labels come out twice and records sit
  # under other records' names. That is what happened, and it looked like
  # records had gone missing.
  #
  # So the name a record carries is left alone and the split is on its run.
  groups <- !all(is.na(frame$group))
  frame$run <- aev_runs(frame$group)
  run_names <- frame$group[!duplicated(as.integer(frame$run))]

  label_lines <- aev_label_lines(frame$name)

  drawable <- frame[frame$status == "ok", , drop = FALSE]

  # A DIFFERENT SET OF RECORDS FROM THE PANEL ABOVE, deliberately. A cell with
  # no deaths has no position on a log scale but a perfectly good residual, so
  # it appears below and not above. The panels are not two views of one list.
  residual_frame <- frame[is.finite(frame$deviance_residual), , drop = FALSE]
  residual_frame$panel <- as_panel(aev_panels[2], nrow(residual_frame))

  # Drawn rather than themed, so they stop at the A/E panel. See `theme_aev()`.
  grid_lines <- data.frame(
    panel = as_panel(aev_panels[1], length(log_breaks)),
    y = log_breaks
  )

  layers <- list(
    # NOTHING IS DRAWN FROM THESE TWO, they only say how far each panel reaches.
    #
    # The first gives every record its x position in every panel, so a record
    # that cannot be drawn still takes up its slot -- five records of which one
    # was drawable once came out as a chart of one. It replaces a coordinate
    # limit, which could not be used once the groups became a second facet
    # dimension, because a coordinate limit is one limit for the whole plot and
    # each group column needs its own.
    #
    # The second gives each panel a y range even when it holds no drawable
    # records at all, so that `aev_panel_limits()` is never handed a degenerate
    # range and left guessing which panel it is looking at.
    ggplot2::geom_blank(
      data = aev_skeleton(frame, panels),
      mapping = ggplot2::aes(x = .data$position)
    ),
    ggplot2::geom_blank(
      data = aev_reach(panels, log_range),
      mapping = ggplot2::aes(y = .data$y)
    ),
    # THE GROUP MUST ARRIVE AS THE FACTOR, not as its label. A character value
    # here makes the combined facet variable a character, and the columns come
    # out in alphabetical order instead of record order -- which is precisely
    # what happened, and what `test-aev-facet.R` caught.
    aev_legend_layer(
      residuals,
      panel = as_panel(aev_panels[1]),
      run = if (nrow(frame) == 0L) NA else frame$run[1]
    ),
    ggplot2::scale_shape_manual(
      values = seq_along(aev_legend_entries(residuals)),
      name = NULL
    ),
    ggplot2::geom_hline(
      data = grid_lines,
      mapping = ggplot2::aes(yintercept = .data$y),
      colour = aev_palette$grid_line,
      linewidth = line_width(aev_geometry$grid_width_px)
    ),
    # THE ONE LINE EVERY MARKER IS READ AGAINST, so it is drawn thick and in
    # colour rather than as a twelfth grid line. It goes down early, behind
    # everything, and shows through the transparent middle of a marker.
    ggplot2::geom_hline(
      data = data.frame(panel = as_panel(aev_panels[1]), y = 0),
      mapping = ggplot2::aes(yintercept = .data$y),
      colour = aev_palette$reference_line,
      linewidth = line_width(aev_geometry$reference_width_px)
    ),
    geom_aev_interval(
      data = drawable,
      mapping = ggplot2::aes(
        x = .data$position,
        y = .data$clamped,
        centre = .data$log_ratio,
        sigma = .data$log_A_over_E_stddev
      )
    ),
    geom_aev_chevron(
      data = aev_off_scale(frame, log_range),
      mapping = ggplot2::aes(x = .data$position, top = .data$top)
    )
  )

  if (residuals) {
    layers <- c(
      list(
        ggplot2::geom_rect(
          data = aev_residual_bands(as_panel(aev_panels[2])),
          mapping = ggplot2::aes(
            ymin = .data$lower,
            ymax = .data$upper,
            alpha = .data$alpha
          ),
          xmin = -Inf,
          xmax = Inf,
          fill = aev_palette$residual_band
        ),
        ggplot2::scale_alpha_identity()
      ),
      layers,
      list(
        ggplot2::geom_hline(
          data = data.frame(panel = as_panel(aev_panels[2]), y = 0),
          mapping = ggplot2::aes(yintercept = .data$y),
          colour = aev_palette$axis_line,
          linewidth = line_width(aev_geometry$border_width_px)
        ),
        geom_aev_residual(
          data = residual_frame,
          mapping = ggplot2::aes(x = .data$position, y = .data$deviance_residual)
        )
      )
    )
  }

  ggplot2::ggplot() +
    layers +
    facet_aev(groups, label_lines, run_names) +
    ggplot2::scale_x_continuous(
      breaks = frame$position,
      # THE NEWLINES ARE REAL, and they are here so that ggplot2 measures the
      # depth of the axis correctly. Where the breaks fall is decided again when
      # the chart is drawn and the column width is finally known. See *Wrapping
      # the record names*.
      labels = aev_wrap_labels(frame$name, label_lines),
      # Half a record's width at each end, which is slide 2's rule for the gap
      # at the edges of the plot area.
      expand = ggplot2::expansion(add = 0.5)
    ) +
    # THE Y LIMITS CANNOT LIVE ON THE COORDINATE SYSTEM, because a coordinate
    # limit is one limit for the whole plot and the two panels need different
    # ones. They go on the scale instead, as a function of each panel's own
    # data range -- which is what `scales = "free_y"` makes possible.
    #
    # `oob` IS WHAT KEEPS THAT SAFE. A scale limit censors by default, turning
    # anything outside into NA and taking whole shapes with it: an interval
    # whose cap runs off the top would vanish entirely rather than be trimmed.
    # Passing values through unchanged leaves the clipping to the panel, which
    # is what the coordinate limit used to do.
    ggplot2::scale_y_continuous(
      limits = function(range) aev_panel_limits(range, log_range),
      oob = function(values, range) values,
      breaks = function(limits) aev_panel_breaks(limits, log_breaks),
      labels = aev_panel_labels,
      expand = ggplot2::expansion(0, 0)
    ) +
    ggplot2::labs(
      title = title, subtitle = subtitle, x = NULL, y = NULL,
      # THE NOTES ARE THE CAPTION. Setting it here rather than leaving it empty
      # means `+ ggplot2::labs(caption = ...)` still wins afterwards, which is
      # the behaviour a ggplot2 user will expect of anything they add with `+`.
      caption = aev_notes_text(notes)
    ) +
    theme_aev() +
    # The residual strip is a fixed depth and the A/E panel takes what is left,
    # which is how `AEChart.cs` lays the two out.
    # The panel heights are NOT set here. See `aev_set_panel_heights()`.
  ggplot2::theme()
}

# One row per record per panel, carrying nothing but where the record sits and
# which group it belongs to.
aev_skeleton <- function(frame, panels) {

  if (nrow(frame) == 0L) return(frame[, c("position", "run", "panel"), drop = FALSE])

  do.call(rbind, lapply(panels, function(name) {
    skeleton <- frame[, c("position", "run"), drop = FALSE]
    skeleton$panel <- factor(rep(name, nrow(frame)), levels = panels)
    skeleton
  }))
}

# Two rows per panel: the ends of the range that panel is drawn at.
aev_reach <- function(panels, log_range = aev_log_ratio_limits) {

  reaches <- list(
    log_range,
    c(-aev_residual$limit, aev_residual$limit)
  )[seq_along(panels)]

  do.call(rbind, Map(function(name, reach) {
    data.frame(panel = factor(rep(name, 2L), levels = panels), y = reach)
  }, panels, reaches))
}

# `unit.c` of a null and an absolute length: the first panel absorbs the slack,
# the second is 84 px at the reference type size and grows with it from there.
aev_panel_heights <- function(residuals, scale = 1) {
  if (!residuals) return(grid::unit(1, "null"))
  grid::unit.c(grid::unit(1, "null"), inches(aev_residual$height_px * scale))
}

#' @rdname aev_plot
#' @export
plot.aev <- function(x, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
