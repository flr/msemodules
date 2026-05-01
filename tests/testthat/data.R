# test-performance-db.R - DESC
# /home/mosqu003/Active/mse_FLR/msemodules/tests/testthat/test-performance-db.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

library(mse)

# BUILD example FLmse and FLmses

data(plesim)

# Set control: sa and hcr
control <- m
library(mse)pCtrl(list(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=hockeystick.hcr, args=list(lim=0,
  trigger=14000, target=0.18))))

# RUN mp
te <- mps(om, oem=oem, ctrl=control, args=list(iy=2021, fy=2034))

# RUN mps
tes <- mps(om, oem=oem, ctrl=control, args=list(iy=2021, fy=2034),
  hcr=list(target=c(0.18, 0.25, 0.25)))

save(te, tes, file="mse.rda", compress="xz")

