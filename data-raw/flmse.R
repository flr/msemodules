# flmse.R - DESC
# /home/mosqu003/Active/mse_FLR/msemodules/data-raw/flmse.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

data(plesim)

# SC devs for SSB
scdevs <- shortcut_devs(om, Fcv=0.212, Fphi=0.423, SSBcv=0.10)

# SETUP control shortcut.sa + hoceystick.hcr + bank_borrow.is

ctrl <- mpCtrl(est = mseCtrl(method=shortcut.sa,
    args=list(SSBdevs=scdevs$SSB)),
  hcr = mseCtrl(method=hockeystick.hcr, args=list(metric="ssb", trigger=9e5, 
    output="catch", target=1.2e5)),
  isys = mseCtrl(method=bank_borrow.is, args=list(rate=0.10, healthy=3, diff=0.01)))

# Runs the MP
run <- mp(om, control=ctrl, args=list(iy=2021, fy=2035, frq=1))

save(run,  file="../data/mserun.rda")

