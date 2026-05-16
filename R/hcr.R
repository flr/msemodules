# hcr.R - DESC
# msemodules/R/hcr.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
#         Ernesto Jardim (IPMA) <ernesto.jardim@ipma.pt>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


globalVariables(c(".", "ay", "bufflow", "buffupp", "data_lag", "data",
  "dy", "frq", "fy", "iy", "lim", "management_lag", "min", "mys",
  "xmin", "xmax", "label", "element", "label", "run", "statistic",
  "sloperatio", "year"))

# buffer.hcr {{{

#' Buffer-based harvest control rule
#'
#' Implements a buffer harvest control rule (HCR) that adjusts management output
#' according to the recent value of a stock or index metric relative to a set of
#' threshold reference points.
#'
#' @param stk A stock object, typically an `FLStock`, used together with `ind`
#'   to compute the HCR metric.
#' @param ind An index object used when computing the selected metric.
#' @param metric Character string specifying the metric to use in the HCR.
#' @param target Numeric. Target value of the metric.
#' @param width Numeric. Half-width of the buffer around `target`.
#' @param lim Numeric. Lower limit reference point for the metric.
#' @param bufflow Numeric. Lower bound of the buffer zone.
#' @param buffupp Numeric. Upper bound of the buffer zone.
#' @param sloperatio Numeric. Controls the slope of output increases above
#'   `buffupp`.
#' @param initial Initial output value used when no previous output is available.
#' @param nyears Number of recent years over which the metric is averaged.
#' @param dupp Optional upper bound on proportional increase in output.
#' @param dlow Optional lower bound on proportional decrease in output.
#' @param all Logical. If `TRUE`, change limits are always applied.
#' @param scale Logical. If `TRUE`, do not recover previous output from tracking.
#' @param ... Additional arguments passed to [selectMetric()].
#' @param args A list of auxiliary arguments unpacked by [spread()].
#' @param tracking A tracking object updated with metric, decision, tier, and
#'   output values.
#' @param x An object for which [args()] returns a list containing at least the
#'   buffer HCR parameters `lim`, `bufflow`, and `buffupp`.
#'
#' @return
#' For [buffer.hcr()], a list with components `ctrl` and `tracking`.
#'
#' For [buffer_bands()], a list of ggplot2 layers for shading HCR zones.
#'
#' @details
#' [buffer.hcr()] computes a piecewise decision multiplier from the selected
#' metric and applies it to the previous output value.
#'
#' [buffer_bands()] creates shaded ggplot2 rectangles marking the critical,
#' recovery, buffer, and above-target zones implied by the same thresholds.
#'
#' @seealso [plot_buffer.hcr()]
#' @name buffer.hcr
#' @author Iago Mosqueira, WMR; Richard Hillary, CSIRO

buffer.hcr <- function(stk, ind, metric='wmean',
  target=1, width=0.5, lim=max(target * 0.10, target - 2 * width), 
  bufflow=max(lim * 1.50, target - width), buffupp=target + width,
  sloperatio=0.15, initial, nyears=4, dupp=NULL, dlow=NULL, all=TRUE, scale=FALSE,
  ..., args, tracking) {

  # EXTRACT args: dy, ay, mys, frq)
  spread(args)

  # COMPUTE and window metric
  met <- yearMeans(window(selectMetric(metric, stk, ind, ...),
    start=dy - nyears + 1, end=dy))

  track(tracking, "metric.hcr", ay) <- met

  # COMPUTE HCR multiplier if ...
  # BELOW lim
  dec <- ifelse(met <= lim, ((met / lim) ^ 2) / 2,
    # BETWEEN lim and bufflow
    ifelse(met <= bufflow,
      (0.5 * (1 + (met - lim) / (bufflow - lim))),
    # BETWEEN bufflow and buffupp
    ifelse(met < buffupp, 1, 
    # ABOVE buffupp, as proportion of downward gradient
      1 + sloperatio * ((1 - 2^(-1)) / (bufflow - lim)) * (met - buffupp))))

  # TRACK rule decision
  track(tracking, "decision.hcr", ay) <- as.numeric(dec)

  # TRACK rule classification (tier)
  tier <- as.numeric(cut(met, c(0, lim, bufflow, buffupp, Inf), labels=seq(1,4)))
  tier[is.na(tier)] <- 0
  track(tracking, "tier.hcr", ay) <- tier

  # GET previous output value if change limited
  
  pre <- FLQuant(initial, iter=args$it)
  
  if(ay > iy & !scale) {
      pre <- tracking[metric == 'hcr' & year == ay - frq, data]
  } 

  # SET TAC
  out <- pre * dec
  
  # TRACK first decision
  track(tracking, "output.hcr", ay) <- out

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
      out[out > pre * dupp] <-
        pre[out > pre * dupp] * dupp
    } else {
      out[out > pre * dupp & met > bufflow] <-
        pre[out > pre * dupp & met > bufflow] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
      out[out < pre * dlow] <-
        pre[out < pre * dlow] * dlow
    } else {
      out[out < pre * dlow & met > bufflow] <-
        pre[out < pre * dlow & met > bufflow] * dlow
    }
  }

  # SUBSTITUTE NAs for 0
  out <- ifelse(is.na(out), 0, out)

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for frq years
    c(lapply(mys, function(x) list(quant="catch", value=c(out), year=x))))
	
  list(ctrl=ctrl, tracking=tracking)
}
# }}}

# plot_buffer.hcr {{{

#' Plot and annotate a buffer harvest control rule
#'
#' Tools for visualizing a buffer-based harvest control rule (HCR), including
#' the HCR curve itself and optional shaded management-tier bands.
#'
#' [plot_buffer.hcr()] draws the HCR curve implied by the buffer rule and can
#' optionally overlay observed values.
#'
#' [buffer_bands()] returns ggplot2 layers that shade the management tiers
#' defined by the HCR thresholds.
#'
#' @param args A list of arguments defining the buffer HCR. Typically includes
#'   values such as `lim`, `bufflow`, `buffupp`, `sloperatio`, and optionally
#'   other elements used by the plotting method. If an object of class
#'   `mseCtrl` is supplied to [plot_buffer.hcr()], its arguments are extracted
#'   using [args()].
#' @param obs Optional observed values to overlay on the HCR plot. Can be:
#'   \itemize{
#'     \item an `FLStock` object, in which case `metric` and `output` are
#'       extracted and plotted;
#'     \item a numeric value, in which case a single point is added on the HCR
#'       curve;
#'     \item omitted, in which case no observations are added.
#'   }
#' @param alpha Numeric. Alpha transparency used for observed points and paths.
#' @param labels Labels for plot annotations. A named vector or list with names
#'   among `lim`, `bufflow`, `buffupp`, `metric`, and `output`. Supplied values
#'   override the defaults.
#' @param metric Character string giving the metric used for the x-axis and,
#'   when `obs` is an `FLStock`, the variable extracted from [metrics()].
#' @param output Character string giving the output used for the y-axis and,
#'   when `obs` is an `FLStock`, the variable extracted from [metrics()].
#' @param xlim Numeric. Maximum x-axis value for the HCR plot. By default set to
#'   `buffupp * 1.5`.
#' @param ylim Numeric. Maximum y-axis value for the HCR plot.
#' @param x An object for which [args()] returns a list containing at least the
#'   buffer-HCR parameters `lim`, `bufflow`, and `buffupp`.
#'
#' @return
#' For [plot_buffer.hcr()], a `ggplot2` object.
#'
#' For [buffer_bands()], a list of ggplot2 components:
#' \itemize{
#'   \item a [ggplot2::geom_rect()] layer that draws shaded tier bands;
#'   \item a [ggplot2::scale_fill_manual()] layer that assigns colours to tiers.
#' }
#'
#' @details
#' The buffer harvest control rule is defined by threshold values that partition
#' the x-axis metric into zones:
#' \itemize{
#'   \item below `lim`: critical zone;
#'   \item between `lim` and `bufflow`: recovery zone;
#'   \item between `bufflow` and `buffupp`: buffer zone;
#'   \item above `buffupp`: above-target zone.
#' }
#'
#' [plot_buffer.hcr()] evaluates the implied HCR multiplier across a sequence of
#' metric values and draws the resulting curve, including annotations for the
#' key thresholds.
#'
#' [buffer_bands()] constructs semi-transparent shaded rectangles spanning these
#' threshold intervals to support interpretation of HCR plots.
#'
#' @section plot_buffer.hcr:
#' `plot_buffer.hcr()` computes the HCR response curve from the supplied buffer
#' rule parameters and returns a ggplot object. The plot includes:
#' \itemize{
#'   \item the HCR line;
#'   \item vertical reference markers for `lim`, `bufflow`, and `buffupp`;
#'   \item labels for the threshold points;
#'   \item optional observed points or trajectories.
#' }
#'
#' If `obs` is an `FLStock` object, the function extracts the chosen `metric`
#' and `output` values using [metrics()] and overlays them on the HCR curve. If
#' there is a single iteration, both points and a path are added, and first/last
#' years are labelled.
#'
#' If `obs` is numeric, a single red point is placed at the corresponding
#' location on the HCR curve.
#'
#' @section buffer_bands:
#' `buffer_bands()` extracts `lim`, `bufflow`, and `buffupp` from `args(x)` and
#' defines four contiguous x-axis intervals. A fourth threshold,
#' `target = 1.5 * buffupp`, is used as the upper bound for the final tier.
#'
#' The returned rectangles extend from `y = 0` to `y = Inf`, so they occupy the
#' full vertical plotting range.
#'
#' @section Assumptions and caveats:
#' \itemize{
#'   \item Thresholds are assumed to satisfy
#'     `0 <= lim <= bufflow <= buffupp`.
#'   \item `buffer_bands()` adds a fill scale and may therefore override an
#'     existing fill scale in a ggplot object.
#'   \item In [plot_buffer.hcr()], default values such as `xlim` and `ylim`
#'     depend on objects created after argument expansion and may rely on the
#'     surrounding evaluation framework.
#'   \item `plot_buffer.hcr()` currently uses internal helper functions such as
#'     [spread()], [args()], and [metrics()], and assumes these are available.
#' }
#'
#' @seealso [args()], [metrics()], [ggplot2::geom_line()],
#'   [ggplot2::geom_rect()], [ggplot2::scale_fill_manual()]
#'
#' @examples
#' \dontrun{
#' args <- list(lim = 0.4, bufflow = 1, buffupp = 2, sloperatio = 0.2)
#'
#' # Plot the HCR curve
#' plot_buffer.hcr(args, labels = list(metric = "CPUE", output = "C~mult"))
#'
#' # Add buffer bands to a custom ggplot
#' layers <- buffer_bands(hcr)
#' p <- ggplot(dat, aes(x = biomass, y = harvest)) + geom_line()
#' for (ly in layers) p <- p + ly
#' p
#' }
#'
#' @author Iago Mosqueira, WMR; Richard Hillary, CSIRO
#' @keywords hplot
#' @name buffer.hcr

plot_buffer.hcr <- function(args, obs="missing", alpha=0.3,
  labels=c(lim="Limit", bufflow="Lower~buffer", buffupp="Upper~buffer",
    metric=metric, output=output), metric=args$metric, output='multiplier',
    xlim=buffupp * 1.50, ylim=scale * 1.50) {

  # EXTRACT args from mpCtrl
  if(is(args, "mseCtrl"))
    args <- args(args)

  # GET args
  spread(lapply(args, c), FORCE=TRUE)

  # PARSE labels
  alllabels <- formals()$labels
  alllabels[names(labels)] <- labels
  labels <- as.list(alllabels)

  # SET met values
  met <- seq(0, xlim, length=200)

  # COMPUTE gradient of decrease
  dgradient <- (1 - 2^(-1))/(bufflow - lim)

  # COMPUTE HCR multiplier if ...

  # BELOW lim
  out <- ifelse(met <= lim, ((met / lim) ^ 2) / 2,
    # BETWEEN lim and bufflow
    ifelse(met <= bufflow,
      (0.5 * (1 + (met - lim) / (bufflow - lim))),
    # BETWEEN bufflow and buffupp
    ifelse(met < buffupp, 1, 
    # ABOVE buffupp, as proportion of dgradient
      1 + sloperatio * dgradient * (met - buffupp))))

  # DATA
  # TODO: ADD 'set'
  dat <- data.frame(met=met, out=out)
  
  # TODO:
  scale <- 1
  
  # TODO: ADD aes(group='set')
  p <- ggplot(dat, aes(x=met, y=out)) +
    coord_cartesian(ylim = c(0, ylim), clip="off") +
    # HCR line
    geom_line() +
    # TODO: TARGET & WIDTH
    # BUFFER UPP
    annotate("segment", x=buffupp, xend=buffupp, y=0, yend=scale, linetype=2) +
    annotate("point", x=buffupp, y=scale, size=3) +
    annotate("text", x=buffupp, y=-ylim / 40, label=labels$buffupp,
      vjust="bottom", parse=TRUE) +
    # BUFFER LOW
    annotate("segment", x=bufflow, xend=bufflow, y=0, yend=scale, linetype=2) +
    annotate("point", x=bufflow, y=scale, size=3) +
    annotate("text", x=bufflow, y=-ylim / 40, label=labels$bufflow,
      vjust="bottom", parse=TRUE) +
    # LIMIT
    annotate("segment", x=lim, xend=lim, y=0, yend=out[which.min(abs(met - lim))], 
      linetype=2) +
    annotate("point", x=lim, y=out[which.min(abs(met - lim))], size=3) +
    annotate("text", x=lim, y=-ylim / 40, label=labels$lim, vjust="bottom",
      parse=TRUE) +
    # SLOPE
    annotate("segment", x=buffupp, xend=xlim, y=1, yend=1, linetype=2) +
    annotate("text", x=buffupp + (xlim - buffupp) / 3, y=1, label="slope", 
      vjust="bottom", parse=TRUE)

  # AXIS labels
  if(!is.null(labels$metric))
    p <- p + xlab(parse(text=labels$metric))
  if(!is.null(labels$output))
    p <- p + ylab(parse(text=labels$output))

  # OBS
  if(!missing(obs)) {
    # FLStock
    if(is.FLStock(obs)) {
      obs <- model.frame(metrics(obs, list(met=get(metric), out=get(output))))
      xlim <- max(obs$met, na.rm=TRUE) * 1.05
      ylim <- max(obs$out, na.rm=TRUE) * 1.05

      # PLOT line if 1 iter
      if(length(unique(obs$iter)) == 1)
        p <- p + geom_point(data=obs, alpha=alpha) +
          geom_path(data=obs, alpha=alpha) +
          geom_label(data=subset(obs, year %in% c(min(year), max(year))),
            aes(label=year), fill=c('gray', 'white'), alpha=1)
      # PLOT with alpha if multiple
      else
        p <- p + geom_point(data=obs, alpha=alpha)
    }
    # NUMERIC
    else if(is.numeric(obs)) {
      obs <- data.frame(met=obs, out=out[which.min(abs(met - obs))])
      p <- p + geom_point(data=obs, colour="red", size=3)
    }

  }
  return(p)
}

# }}}

# buffer_bands {{{

#' @rdname buffer.hcr

buffer_bands <- function(x) {

  # ARGS
  lim <- args(x)$lim
  buffupp <- args(x)$buffupp
  bufflow <- args(x)$bufflow
  target <- buffupp * 1.5

  tier_bands <- data.frame(
    xmin  = c(0,       lim,     bufflow, buffupp),
    xmax  = c(lim,     bufflow, buffupp, target),
    label = c("Tier 1\n(critical)", "Tier 2\n(recovery)",
            "Tier 3\n(buffer)",   "Tier 4\n(above target)"))

  tier_colours <- c(
    "Tier 1\n(critical)"      = "#A32D2D",
    "Tier 2\n(recovery)"      = "#854F0B",
    "Tier 3\n(buffer)"        = "#3B6D11",
    "Tier 4\n(above target)"  = "#185FA5")

  return(list(geom_rect(data = tier_bands,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = label),
      inherit.aes = FALSE, alpha = 0.4),
    scale_fill_manual(values = tier_colours, name = NULL)))
}

# }}}

# depletion.hcr {{{

#' Depletion-based harvest control rule
#'
#' Implements a depletion-based harvest control rule (HCR) in which management
#' output is scaled according to stock status relative to carrying capacity.
#'
#' The rule uses the ratio of a selected stock metric to `K` as a depletion
#' indicator and adjusts the harvest rate relative to `hrmsy` depending on
#' whether the stock is above a trigger level, between a limit and trigger, or
#' below the limit.
#'
#' @param stk A stock object, typically an `FLStock`, used together with `ind`
#'   to compute the selected metric.
#' @param ind An FLQuant object containing metrics to be used by the HCR.
#' @param metric Character string specifying the metric to use in the HCR.
#'   Passed to [selectMetric()]. Defaults to `"ssb"`.
#' @param mult Numeric multiplier applied to the harvest rate after scaling by
#'   the decision rule. Defaults to `1`.
#' @param hrmsy Numeric target harvest rate associated with MSY conditions.
#' @param K Numeric carrying capacity or proxy scaling constant used to convert
#'   the metric into a depletion ratio.
#' @param trigger Numeric depletion threshold above which the HCR applies the
#'   full harvest rate. Defaults to `0.4`.
#' @param lim Numeric depletion limit reference point below which the HCR
#'   applies the minimum multiplier. Defaults to `0.1`.
#' @param min Numeric minimum multiplier applied when depletion falls below
#'   `lim`. Defaults to `0.00001`.
#' @param initial Optional numeric initial output value used when interannual
#'   change constraints are applied in the first management year.
#' @param dupp Optional numeric upper bound on proportional increase in output.
#'   If not `NULL`, output increases are capped at `pre * dupp`.
#' @param dlow Optional numeric lower bound on proportional decrease in output.
#'   If not `NULL`, output decreases are floored at `pre * dlow`.
#' @param all Logical. If `TRUE`, change constraints are always applied. If
#'   `FALSE`, they are only applied when the metric is below `trigger`.
#'   Defaults to `TRUE`.
#' @param ... Additional arguments intended for downstream metric calculation.
#' @param args A list of auxiliary framework arguments unpacked by [spread()].
#'   The function expects at least `ay`, `iy`, `dy`, `mys`, `data_lag`,
#'   `management_lag`, and `it`.
#' @param tracking A tracking object updated throughout the [mse::mp()] evaluation.
#'
#' @return
#' A named list with components:
#' \itemize{
#'   \item `ctrl`: an `fwdControl` object specifying the catch target for the
#'     management years in `mys`;
#'   \item `tracking`: the updated tracking object containing recorded metric,
#'     decision, tier, and output values.
#' }
#'
#' @details
#' The function first computes the selected metric from `stk` and `ind` using
#' [selectMetric()], restricted to the terminal data year `dy`. This metric is
#' tracked under `"metric.hcr"`.
#'
#' Stock status is then expressed as the depletion ratio:
#' \deqn{met / K}
#'
#' The decision multiplier is calculated as a piecewise function:
#' \itemize{
#'   \item if `met / K >= trigger`, the multiplier is `1`;
#'   \item if `lim <= met / K < trigger`, the multiplier increases linearly from
#'     `0` at `lim` to `1` at `trigger`;
#'   \item if `met / K < lim`, the multiplier is set to `min`.
#' }
#'
#' This multiplier is recorded in `tracking` under `"decision.hcr"`.
#'
#' The rule also classifies stock status into three tiers based on depletion:
#' \itemize{
#'   \item Tier 1: below `lim`;
#'   \item Tier 2: between `lim` and `trigger`;
#'   \item Tier 3: at or above `trigger`.
#' }
#'
#' Tier membership is stored under `"tier.hcr"`.
#'
#' A harvest-rate target is then computed as:
#' \deqn{hrtarget = hrmsy \times dec \times mult}
#'
#' and management output is set as:
#' \deqn{out = met \times hrtarget}
#'
#' This provisional output is stored under `"output.hcr"` before optional upper
#' and lower change constraints are applied.
#'
#' If `dupp` and/or `dlow` are supplied, the previous output value is recovered
#' either from `initial` in the first management year or from the tracking
#' object in later years. The provisional output is then capped and/or floored
#' according to the specified proportional limits.
#'
#' Finally, missing values are replaced with zero and the result is converted to
#' an [fwdControl()] object with one catch target per management year.
#'
#' @section Change constraints:
#' If `dupp` is supplied, increases in output can be capped. If `dlow` is
#' supplied, decreases can be limited. When `all = TRUE`, these constraints are
#' always applied. When `all = FALSE`, they are only applied when the metric is
#' below the `trigger`, allowing unconstrained changes when stock status is at
#' or above the trigger level.
#'
#' @section Assumptions and caveats:
#' \itemize{
#'   \item The function assumes `spread(args[...])` creates `ay`, `iy`, `dy`,
#'     `mys`, and `it` in scope.
#'   \item Although `...` is declared, it is not currently passed on in the call
#'     to [selectMetric()] in the implementation shown.
#'   \item The previous output retrieval uses
#'     `tracking[metric == 'hcr' & year == ay, data]`, which assumes a specific
#'     internal structure for `tracking`.
#'   \item `track(tracking, "decision.hcr", mys)` stores the decision for all
#'     management years, whereas `track(tracking, "tier.hcr", ay)` stores tier
#'     at assessment year resolution.
#'   \item Missing outputs are silently converted to zero.
#' }
#'
#' @seealso [buffer.hcr()], [selectMetric()], [fwdControl()], [mse::`track()<-`]
#'
#' @author Iago Mosqueira, WMR
#' @keywords models

depletion.hcr <- function(stk, ind, metric='ssb', mult=1, hrmsy, K,
  trigger=0.4, lim=0.1, min=0.00001, initial, dupp=NULL, dlow=NULL, all=TRUE,
  ..., args, tracking) {

  # EXTRACT args
  spread(args[c('ay', 'iy', 'dy', 'mys', 'data_lag', 'management_lag', 'it')])

  # COMPUTE and window metric
  met <- window(selectMetric(metric, stk, ind), start=dy, end=dy)

  # TRACK metric
  track(tracking, "metric.hcr", dy) <- met

  # COMPUTE HCR multiplier if ...
  dec <- ifelse((met / K) >= trigger, 1,
    # metric betwween lim and trigger, and above
    ifelse((met / K) >= lim,      
      ((met / K) - lim) / (trigger - lim), min))

  # TRACK initial output
  track(tracking, "decision.hcr", mys) <- dec

  # TRACK decision
  tier <- cut((met / K), c(0, lim, trigger, Inf), labels=seq(1, 3))

  track(tracking, "tier.hcr", ay) <- as.numeric(tier)

  # SET HR target
  hrtarget <- hrmsy * dec * mult

  # SET TAC
  out <- met * hrtarget

  # TRACK first decision
  track(tracking, "output.hcr", mys) <- out

  # GET previous output value if change limited
  if(!is.null(dupp) | !is.null(dlow)) {
    # GET initial value at start if set,
    if(ay == iy) {
      # STOP if initial is NULL
      if(is.null(initial))
        stop("To apply 'dlow' and 'dupp' limits, 'initial' is required")
      pre <- FLQuant(initial, iter=args$it)
    # OR previous decision,
    } else {
      pre <- tracking[metric == 'hcr' & year == ay, data]
    }
  }

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
      out[out > pre * dupp] <-
        pre[out > pre * dupp] * dupp
    } else {
      out[out > pre * dupp & met < trigger] <-
        pre[out > pre * dupp & met < trigger] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
      out[out < pre * dlow] <-
        pre[out < pre * dlow] * dlow
    } else {
      out[out < pre * dlow & met < trigger] <-
        pre[out < pre * dlow & met < trigger] * dlow
    }
  }

  # SUBSTITUTE NAs for 0
  out <- ifelse(is.na(out), 0, out)

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for frq years
    c(lapply(mys, function(x) list(quant="catch", value=c(out), year=x))))
	
  list(ctrl=ctrl, tracking=tracking)
}
# }}}

# pid.hcr {{{

#' Proportional-Integral-Derivative (PID) Harvest Control Rule (HCR)
#'
#' Implements a Proportional-Integral-Derivative (PID) based harvest control rule for adjusting 
#' Total Allowable Catch (TAC) based on divergence from a reference point using PID control signals.
#'
#' @param stk `FLStock`. The stock object to which the HCR applies.
#' @param ind `FLQuant`. The abundance index used to compute the control signal.
#' @param ref Numeric. The reference value for the metric to determine divergence.
#' @param metric Character or function. The metric applied to measure stock status. Default is `'ssb'` (spawning stock biomass).

pid.hcr <- function(stk, ind, ref, metric=ssb, initial, kp=0, ki=0, kd=0,
  nyears=5, dlow=NA, dupp=NA, args, tracking, ...) {
  
  # args
  spread(args)
  
  # SELECT metric
  met <- window(selectMetric(metric, stk, ind, ...), start=dy - nyears + 1, end=dy)

  # TRACK metric
  track(tracking, "metric.hcr", dy) <- met[, ac(dy)]
 
  # GET TAC from tracking['hcr',]
  if(ay > iy)
    initial <- tracking[[1]]['hcr', ac(ay)]

  # CALCULATE divergence
  # e <- log(met[,5] %/% met[, 4])
  # e <- pmin(pmax(log(met[,2:5] %/% met[, 1:4]), 0.8), 1.2)
  e <- log(met %/% FLQuant(ref))

  # COMPUTE control signal
  u <- kp * e[, ac(dy)] + ki * yearSums(e) + kd * (e[, ac(dy)] - e[, ac(dy - 1)])

  track(tracking, "decision.hcr", ay) <- u

  # COMPUTE factor
  fac <- min(max(exp(u^1), exp(u^2)), exp(u^3))

  # TAC, set to 0 if NA
  tac <- ifelse(is.na(fac * initial), 0, fac * initial)

  # TRACK initial TAC
  track(tracking, "tac.hcr", mys) <- tac

  # LIMITS over previous output
  if(!is.na(dupp))
    tac[tac > initial * dupp] <- initial[tac > initial * dupp] * dupp
  if(!is.na(dlow))
    tac[tac < initial * dlow] <- initial[tac < initial * dlow] * dlow

  # CONTROL
  ctrl <- fwdControl(
    # TAC for frq years
    lapply(seq(ay + management_lag, ay + frq), function(x)
    list(quant="catch", value=c(tac), year=x))
  )

  return(list(ctrl=ctrl, tracking=tracking))
}

# }}}

# cpue.hcr {{{

cpue.hcr <- function(stk, ind, k1=0.2, k2=0.2, k3=0.2, k4=0.2, target=1,
  dtaclow=0.85, dtacupp=1.15, initac=NULL, slope="slope", mean="mean",
  args, tracking) {

  # args
  ay <- args$ay
  frq <- args$frq
  man_lag <- args$management_lag

  # RECOVER slope & mean(cpue)
  slope <- ind[[slope]]
  mcpue <- ind[[mean]]

  # CALCULATE new tac
  ka <- ifelse(slope > 0, k1, k2)
  kb <- ifelse(mcpue > target, k3, k4)

  # GET previous TAC from last hcr ...
  if(is.null(initac)) {
    tac <- areaSums(seasonSums(tracking[[1]]['hcr', ac(ay)]))
    # ... OR catch
    if(all(is.na(tac)))
      tac <- areaSums(seasonSums(catch(stk)[, ac(ay - args$data_lag)]))
  } else {
    tac <- initac
  }

  # TAC_y-1 ~ TAC_y * 1 + ka * m + kb * (mcpue - target)
  tac <- tac * (1 + ka * slope + kb * (mcpue - target))
  
  # TODO: TAC limits, not on 1st year(s)

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for frq years
    c(lapply(seq(ay + man_lag, ay + frq), function(x)
      list(quant="catch", value=c(tac), year=x)))
  )

	return(list(ctrl=ctrl, tracking=tracking))
} # }}}
