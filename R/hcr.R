# hcr.R - DESC
# msemodules/R/hcr.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
#         Ernesto Jardim (IPMA) <ernesto.jardim@ipma.pt>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


globalVariables(c(".", "ay", "bufflow", "buffupp", "data_lag", "data",
  "dy", "frq", "fy", "iy", "lim", "management_lag", "min", "mys",
  "element", "label", "run", "statistic", "sloperatio", "year"))

# buffer.hcr {{{

#' Buffer Harvest Control Rule (HCR)
#'
#' Implements a "buffer" style harvest control rule (HCR) to set changes in catch levels by using a given metric, such as the value of an index of abundance, or an estimate of biomass. This is compared with some defined buffer and limit reference points, to compute a multiplier that will set future catches from previous levels.
#' with optional bounds on TAC changes.
#'
#' @param stk The 'oem' observation or SA estimation, an FLStock object.
#' @param ind Possible indicators returned by the 'est' step, FLQuants.
#' @param metric Character or function. The metric applied by theb rule as input.Default is `wmean`, returned by `cpue.ind`.
#' @param target The desired target value for the metric, defaults to 1. Numeric or FLQuant.
#' @param width Numeric. The width of the buffer zone surrounding the target. Default is `0.5`.
#' @param lim Numeric. The lower threshold beyond which the HCR steeply reduces catch. Default is `max(target * 0.10, target - 2 * width)`.
#' @param bufflow Numeric. The lower bound of the buffer range. Default is `max(lim * 1.50, target - width)`.
#' @param buffupp Numeric. The upper bound of the buffer range. Default is `target + width`.
#' @param sloperatio Numeric. Defines the slope for values exceeding the buffer's upper limit. Default is `0.15`.
#' @param initial Numeric or `FLQuant`. The previous year's TAC value or another scaling factor for outputs.
#' @param nyears Numeric. Number of years to use for windowing or smoothing metrics. Default is 4.
#' @param dlow A limit for the decrease in the output variable, e.g. 0.85 for a maximum decrease of 15%, numeric. No limit is applied if NULL.
#' @param dupp A limit for the increase in the output variable, e.g. 1.15 for a maximum increase of 15%, numeric. No limit is applied if NULL.
#' @param all If `TRUE`, upper and lower limits (`dupp` and `dlow`) are applied unconditionally, otherwise only when metric > bufflow, logical.
#' @param ... Any extra arguments to be passed to the function computing 'metric'.
#' @param args A list containing dimensionality arguments, passed on by mp().
#' @param tracking A data.table used for tracking indicators, intermediate values, and decisions during MP evaluation.

buffer.hcr <- function(stk, ind, metric='wmean',
  target=1, width=0.5, lim=max(target * 0.10, target - 2 * width), 
  bufflow=max(lim * 1.50, target - width), buffupp=target + width,
  sloperatio=0.15, initial, nyears=4, dupp=NULL, dlow=NULL, all=TRUE, scale=FALSE,
  ..., args, tracking) {

  # EXTRACT args
  spread(args, c("dy", "ay", "mys"))

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


  # GET starting output
  pre <- FLQuant(initial, iter=args$it)

  # GET previous output value if change limited
  if(ay > iy & !scale) {
      pre <- tracking[metric == 'hcr' & year == ay, data]
  }

  # SET TAC
  out <- pre * dec
  
  # TRACK first decision
  track(tracking, "output.hcr", mys) <- out

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

#' @rdname buffer.hcr
#' @details
#' `plot_buffer.hcr` Plots the buffer harvest control rule (HCR) curve along with 
#' optional observed values.
#' @param args list. The list of parameters used by the buffer.hcr function.
#' @param obs  Observed values to overlay on the plot.
#' @param alpha Alpha transparency for observed points/lines.
#' @param labels Labels for the plot annotations. A list with names among 'lim', 'bufflow', 'buffupp', 'metric', and 'output'.
#' @param metric The metric name used for the x-axis label.
#' @param output The output name used for the y-axis label.
#' @param xlim The maximum x-axis limit for the plot.
#' @param ylim The maximum y-axis limit for the plot.
#' @return ggplot2 object

plot_buffer.hcr <- function(args, obs="missing", alpha=0.3,
  labels=c(lim="limit", bufflow="Lower~buffer", buffupp="Upper~buffer",
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

# args <- list(lim=0.4, bufflow=1, buffupp=2,
#   sloperatio=0.2)

# plot_buffer.hcr(args, labels=list(metric='CPUE', output='C~mult'))

# }}}

# depletion.hcr {{{

#' Depletion Harvest Control Rule (HCR)
#'
#' Implements a depletion-based harvest control rule (HCR), adjusting harvest rates and TAC (Total Allowable Catch)
#' based on stock status relative to depletion thresholds.
#'
#' @param stk `FLStock`. The stock object to which the HCR applies.
#' @param ind `FLQuant`. The abundance index used as input for computing adjustments.
#' @param metric Character or function. The metric applied to measure stock status. Default is `'ssb'` (spawning stock biomass).
#' @param mult Numeric. A scaling multiplier for adjusting the harvest rate. Default is `1`.
#' @param hrmsy Numeric. Harvest rate that achieves maximum sustainable yield (MSY

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
