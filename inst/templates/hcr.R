# hcr.R - Template for a 'hcr' module for package mse

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


template.hcr <- function(stk, ind, metric, output, %MODULEARGS%,  args, tracking) {

  # EXTRACT args
  spread(args)

  # COMPUTE metric
  # met <- window(selectMetric(metric, stk, ind), start=dy, end=dy)

  # TRACK metric
  track(tracking, "metric.hcr", dy) <- met
  
  # - APPLY rule

  # BELOW lim
  out <- c(ifelse(met <= lim, min,
    # BETWEEN lim and trigger
    ifelse(met < trigger,
      # diff(met - lim) * gradient + min
      (met - lim) * ((target - min) / (trigger - lim)) + min,
    # ABOVE trigger
    target)))

  # TRACK initial output
  track(tracking, paste0(output, ".hcr"), mys) <- out

  # TRACK decision: met <= lim, 1; lim < met < trigger, 2; met >= trigger, 3
  track(tracking, "decision.hcr", ay) <- ifelse(met < drop, 0,
    ifelse(met <= lim, 1, ifelse(met < trigger, 2, 3)))

  # GET TAC dy / ay - 1
  if(ay == iy) {
    pre <- areaSums(unitSums(seasonSums(window(do.call(output, list(stk)),
      start=ay - management_lag, end=ay - management_lag))))
  } else {
    pre <- c(tracking[[1]]["hcr", ac(ay)])
  }

  # IF NA, set to previous value
  if(any(is.na(out))) {
    out[is.na(out)] <- pre[is.na(out)]
  }

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
    out[out > pre * dupp] <- pre[out > pre * dupp] * dupp
    } else {
    out[out > pre * dupp & met < trigger] <- pre[out > pre * dupp & met <
      trigger] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
    out[out < pre * dlow] <- pre[out < pre * dlow] * dlow
    } else {
    out[out < pre * dlow & met < trigger] <- pre[out < pre * dlow & met <
      trigger] * dlow
    }
  }

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for mys years
    c(lapply(mys, function(x) list(quant=output, value=c(out), year=x)))
  )

  # SET fbar ages
  if(output %in% c("f", "fbar")) {
    ctrl$minAge <- range(stk, "minfbar")
    ctrl$maxAge <- range(stk, "maxfbar")
  }

	list(ctrl=ctrl, tracking=tracking)
}
# }}}

