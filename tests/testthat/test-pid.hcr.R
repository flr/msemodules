# test-pid.hcr.R - DESC
# /home/mosqu003/Projects/FLR/code/mse/msetools/msetools/tests/testthat/test-pid.hcr.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(msetools)

data(sol274)

ctrl <- mpCtrl(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=pid.hcr,
   args=list(metric=ssb, ref=refpts(om)$SBMSY, kp=0.5, ki=0.01, kd=0.7)))

run <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2017, fy=2019), .DEBUG=TRUE)
run <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2017, fy=2026))

plot(om, PID=run)


ctrl <- mpCtrl(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=pid.hcr,
   args=list(metric=depletion, ref=0.40, kp=0, ki=0, kd=0,
      initial=15000, B0=2.4e5)))

run <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2017, fy=2032))

plot(om, run)


ctrl <- mpCtrl(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=pid.hcr,
   args=list(metric='depletion', ref=0.40, kp=0.5, ki=0.01, kd=0.7,
      initial=15000, B0=1e5)))


run <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2017, fy=2032))

# C=0

ctrl <- mpCtrl(est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=pid.hcr,
   args=list(metric=ssb, ref=41000, kp=0, ki=0, kd=0, initial=15000)))

run <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2017, fy=2022))

catch(run) 

plot(om, run) + geom_vline(xintercept=2017)

# Cy = Cy-1

ctrl <- mpCtrl(est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=pid.hcr,
   args=list(metric=ssb, ref=41000, kp=1, ki=0, kd=0, initial=15000)))

run <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2017, fy=2022))

catch(run) 

plot(om, run) + geom_vline(xintercept=2017)





# TUNE

# - SET ki = kd = 0
# - INCREASE Kp until output oscillates w/constant amplitude
# - COMBINE Kp (critical gain) with oscillation period (HOW?)
# - ARRIVE at standard tunings for kp, ki and/or kd
#   - strictly proportional (P)
#   - proportional-integral (PI)
#   - proportional-integral-derivative (PID) control
# - kp: proportional response given the control signal.
# - ki: integral part of the response given the sum of recent control signals.
# - kd: derivative part of the response driven by the rate of change in control signal.

