# ind.R - DESC
# msemodules/R/ind.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# len.ind {{{

# TODO

# - SD vs. age, does it increase?
# sd = len * cv
# TODO: arg for slot
# 1. lenSamples(metric(oem))
# 2. metric(lenSamples(perfect.oem), selex)

#' @examples
#' data(ple4)
#' data(ple4.indices)
#' len.ind(ple4, ple4.indices, args=list(ay=2018, data_lag=1),
#'  tracking=FLQuant(dimnames=list(year=2018, metric='lend.ind')),
#'  params=FLPar(linf=132, k=0.080, t0=-0.35))
#' #
#' len.ind(ple4, ple4.indices, args=list(ay=2018, data_lag=1),
#'  indicators=c(lbar=lbar, lmean),
#'  tracking=FLQuant(dimnames=list(year=2018, metric='len.ind')),
#'  params=FLPar(linf=132, k=0.080, t0=-0.35))

# TODO: LOOP over metric

len.ind <- function (stk, idx, args, tracking, indicators="lbar", params,
  nyears=3, cv=0.1, lmax=1.25, bin=1, n=500,
  metric=function(stk) catch.n(stk), ...) {

  # EXTRACT args
  ay <- args$ay
  data_lag <- args$data_lag
  args0 <- list(...)

  # COMPUTE inverse ALK (cv, lmax, bin)
  ialk <- invALK(params, age=seq(dims(stk)$min, dims(stk)$max),
    cv=cv, lmax=lmax, bin=bin)

  # GENERATE length samples from metric
  input <- do.call(metric, list(stk=stk, idx=idx)[names(formals(metric))])

  samps <- lenSamples(window(input, start=ay - data_lag - nyears + 1,
    end=ay - data_lag), ialk, n=n)

  # CONVERT params
  pars <- as(params, "list")
  
  # OBTAIN names from functions
  nms <- unlist(lapply(indicators, function(x)
      if(is(x, "function"))
        find_original_name(x)
      else
        x
    ))

  # SORT OUT names
  if(is.null(names(indicators)))
    names(indicators) <- nms
  else
    names(indicators)[names(indicators) == 0] <- nms[names(indicators) == 0]

  # COMPUTE indicator(s)
  ind <- lapply(setNames(indicators, nm=nms), function(x) {

    # SUBSET indicator arguments in params
    if(is.character(x)) {
      pars <- pars[names(pars) %in% names(formals(get(x)))]
    } else if(is.function(x)) {
      pars <- pars[names(pars) %in% names(formals(x))]
    }

    return(do.call(x, args=c(list(samps), args0, pars)))
  })

  # ADD to tracking on 'ay' as mean across years
  track(tracking, "len.ind", ac(ay)) <- yearMeans(ind[[1]])

  list(stk = stk, ind = FLQuants(ind), tracking = tracking)
}
# }}}
