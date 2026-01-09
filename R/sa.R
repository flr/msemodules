# sa.R - Status and trend estimators and indicators.
# msemodules/R/sa.R

# Copyright European Union, 2018-2021
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# shortcut.sa {{{

shortcut.sa <- function(stk, idx, metric="ssb", SSBdevs=ind %=% 1, devs=SSBdevs,
  args, tracking, ...) {

  # DIMS
  spread(args, names=c('y0', 'dy', 'ay'))

  # SUBSET oem stock
  stk <- window(stk, end=dy)

  # COMPUTE 'metric'
  met <- do.call(metric, c(list(stk), list(...)))

  # CREATE ind from metric and deviances
  ind <- FLQuants(met * window(devs, start=y0, end=dy))

  # NAME as metric
  names(ind) <- metric

  # TRACK 'convergence'
  track(tracking, "conv.est", ay) <- 1
 
  list(stk=stk, ind=ind, tracking=tracking)
}

# }}}

# shortcut_devs {{{

shortcut_devs <- function(om, Fcv=0.212, Fphi=0.423, SSBcv=0, SSBphi=0,
  bias.correct=FALSE) {

  devs <- FLQuants(
    F=rlnormar1(n=dims(om)$iter, meanlog=0, sdlog=Fcv, rho=Fphi, 
      years=dimnames(om)$year, bias.correct=bias.correct),
    SSB=rlnormar1(n=dims(om)$iter, meanlog=0, sdlog=SSBcv, rho=0,
      years=dimnames(om)$year, bias.correct=bias.correct)
  )

  return(devs)
}
# }}}
