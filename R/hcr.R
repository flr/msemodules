# hcr.R - DESC
# msemodules/R/hcr.R

# Copyright European Union, 2018
# Author: Ernesto Jardim (EC JRC) <ernesto.jardim@ec.europa.eu>
#         Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


globalVariables(c("ay", "bufflow", "buffup", "data_lag", "dy", "frq", "fy", "lim",
  "management_lag", "min", "sloperatio"))

# buffer.hcr {{{

buffer.hcr <- function(stk, ind, metric='wmean',
  target=1, width=0.5, lim=max(target * 0.10, target - 2 * width), 
  bufflow=maax(lim * 1.50, target - width), buffupp=target + width,
  sloperatio=0.15, initial, nyears=4, dupp=NULL, dlow=NULL, all=TRUE,
  ..., args, tracking) {

  # EXTRACT args
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

  # TRACK rule classification
  tier <- as.numeric(cut(met, c(0, lim, bufflow, buffupp, Inf), labels=seq(1,4)))
  tier[is.na(tier)] <- 0
  
  track(tracking, "tier.hcr", ay) <- tier

  # GET previous TAC
  if(ay > iy) {
     initial <- tracking[[1]]['hcr', ac(ay)]
  }

  # SET TAC
  out <- initial * dec
  
  # TRACK first decision
  track(tracking, "output.hcr", mys) <- out

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
      out[out > initial * dupp] <-
        initial[out > initial * dupp] * dupp
    } else {
      out[out > initial * dupp & met > bufflow] <-
        initial[out > initial * dupp & met > bufflow] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
      out[out < initial * dlow] <-
        initial[out < initial * dlow] * dlow
    } else {
      out[out < initial * dlow & met > bufflow] <-
        initial[out < initial * dlow & met > bufflow] * dlow
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

depletion.hcr <- function(stk, ind, metric='ssb', mult=1, hrmsy, K,
  trigger=0.4, lim=0.1, min=0.00001, dupp=NULL, dlow=NULL, all=TRUE,
  ..., args, tracking) {

  # EXTRACT args
  ay <- args$ay
  iy <- args$iy
  data_lag <- args$data_lag
  man_lag <- args$management_lag
  frq <- args$frq

  # SET data year
  dy <- ay - data_lag
  # SET control years
  cys <- seq(ay + man_lag, ay + man_lag + frq - 1)

  # COMPUTE and window metric
  met <- mse::selectMetric(metric, stk, ind)
  met <- window(met, start=dy, end=dy)

  # COMPUTE HCR multiplier if ...
  hcrm <- ifelse((met / K) >= trigger, 1,
    # metric betwween lim and trigger, and above
    ifelse((met / K) >= lim,      
      ((met / K) - lim) / (trigger - lim), min))

  # TRACK decision
  dec <- cut((met / K), c(0, lim, trigger, Inf), labels=seq(1, 3))
  # track(tracking, "decision.hcr", ay) <- as.numeric(dec)

  # SET HR target
  hrtarget <- hrmsy * hcrm * mult

  # SET TAC
  out <- met * hrtarget

  # TRACK first decision
  track(tracking, "rule.hcr", cys) <- out

  # TODO: ADD initac and TAC from tracking
  lastcatch <- unitSums(seasonSums(window(catch(stk), start=dy, end=dy)))

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
      out[out > lastcatch * dupp] <-
        lastcatch[out > lastcatch * dupp] * dupp
    } else {
      out[out > lastcatch * dupp & met < trigger] <-
        lastcatch[out > lastcatch * dupp & met < trigger] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
      out[out < lastcatch * dlow] <-
        lastcatch[out < lastcatch * dlow] * dlow
    } else {
      out[out < lastcatch * dlow & met < trigger] <-
        lastcatch[out < lastcatch * dlow & met < trigger] * dlow
    }
  }

  # TRACK TAC limited
  # sum(old > out)
  # sum(old < out)

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for frq years
    c(lapply(cys, function(x) list(quant="catch", value=c(out), year=x))))
	
  list(ctrl=ctrl, tracking=tracking)
}
# }}}

# pid.hcr {{{

#' @param stk
#' @param ind
#' @param kp
#' @param ki
#' @param kd
#' @param nyears Number of years used in regression of log(stock).
#' @param metric
#' @param ref
#' @examples
#' #control <- mpCtrl(list(
#' #  est = mseCtrl(method=perfect.sa),
#' #  hcr = mseCtrl(method=pid.hcr,
#' #   args=list(metric=ssb, ref=refpts(om)$SBMSY, kp=0.5, ki=0.01, kd=0.7))))
#' # tes <- mp(om, oem=oem, ctrl=control, args=list(iy=2017))
#' # plot(om, PID=tes)

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
    tac[tac > pre * dupp] <- pre[tac > pre * dupp] * dupp
  if(!is.na(dlow))
    tac[tac < pre * dlow] <- pre[tac < pre * dlow] * dlow

  # CONTROL
  ctrl <- fwdControl(
    # TAC for frq years
    lapply(seq(ay + management_lag, ay + frq), function(x)
    list(quant="catch", value=c(tac), year=x))
  )

  return(list(ctrl=ctrl, tracking=tracking))
}

# }}}
