# utilities.R - DESC
# msemodules/R/utilities.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# smooth-index {{{
smooth_index <- function(x, enp.mult=0.2) {

  x[x==0]<-1e3
  
  dat <- data.table(as.data.frame(x, drop=TRUE))
  dat <- dat[!is.na(data)]

  dat[, enptarget:=sum(!is.na(data)) * enp.mult, by=iter]

  dat[, data:=exp(predict(loess(log(data)~year,
    enp.target=unique(enptarget)))), by=iter]

  return(as.FLQuant(dat[, .(year, iter, data, age='all')]))
}
# }}}

