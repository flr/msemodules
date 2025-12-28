# test-buffer.hcr.R - DESC
# msemodules/tests/testthat/test-buffer.hcr.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(msemodules)

data(plesim)

# SET dims

iy <- 2010

# EXAMPLE: perfect.sa + buffer(fbar ~ depletion)

ctrl <- mpCtrl(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=buffer.hcr,
    args=list(metric='depletion', lim=0.20, bufflow=0.35, buffupp=0.55, sloperatio=0.15,
    nyears=1, initial=mean(catch(om)[, ac(iy)]), B0=refpts(om)$B0 * 0.92)))

ctrl <- mpCtrl(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(metric='depletion', lim=0.20, trigger=0.40,
      output="catch", target=mean(catch(om)[, ac(iy)]),
      initial=mean(catch(om)[, ac(iy)]), B0=refpts(om)$B0 * 0.92)))

run <- mp(om, ctrl=ctrl, args=list(iy=iy, frq=1))

runs <- lapply(setNames(nm=c(1990, 2005, 2020)), function(i) {
  mp(om, ctrl=ctrl, args=list(iy=i, frq=1))
})

plot(om, FLmses(runs), window=2020)

run <- mp(iter(om, 1), ctrl=iter(ctrl, 1), args=list(iy=iy, fy=2020, frq=1))

# PLOT run

# PLOT BUG: names, metrics
plot(om, run)

# DEBUG:
# plot(om, run, metrics=list(SB=ssb, C=catch, R=rec, D=depletion))
# metrics(om, metrics=list(SB=ssb, D=depletion))

# PLOT metric (D) & output (C)
plot(FLQuants(D=depletion(stock(run), B0=refpts(om)$B0), C=catch(run)))

# CHECK decisions
# - metric.hcr
# - decision.hcr
# - tier.hcr
# - output.hcr
# - hcr

ays <- args(run)$vy

tra <- FLQuants(
  metric=tracking(run)['metric.hcr', ays],
  decision=tracking(run)['decision.hcr', ays],
  tier=tracking(run)['tier.hcr', ays],
  output=tracking(run)['output.hcr', ays],
  advice=tracking(run)['hcr', ac(an(ays) + 1)],
  status=tracking(run)['SB.om', ac(an(ays) + 1)],
  effect=tracking(run)['C.om', ac(an(ays) + 1)]
)

plot(tra)

tra$output[, ac(1992:2019)] / tra$output[, ac(1991:2018)]
tra$effect[, ac(1992:2019)] / tra$effect[, ac(1991:2018)]


# OBS: met, out
obs <- model.frame(tra[c("metric", "decision")], drop=TRUE)
names(obs) <- c('year', 'met', 'out')

plot_buffer.hcr(args(control(run, 'hcr'))) +
  geom_point(data=obs, colour="red", alpha=0.2)

plot_buffer.hcr(args(control(run, 'hcr'))) +
  geom_path(data=subset(obs, iter==1), colour="red", alpha=0.2) +
  geom_point(data=subset(obs, iter==1), colour="red", alpha=0.6) +
  geom_point(data=subset(obs, iter==1 & year == 1991), size=2)


ext <- FLQuants(
  # ay, from dy
  basis=tracking(run)['metric.hcr', ac(ay)],
  # ay
  rule=tracking(run)['rule.hcr', ac(ay)],
  # ay
  decision=tracking(run)['decision.hcr', ac(ay)],
  # mys
  output=tracking(run)['hcr', ac(ay)],
  # msy
  management=tracking(run)['fwd', ac(ay)],
  # mys
  change=catch(run)[, ac(mys)] / catch(run)[, ac(mys - 1)],
  # ay
  status=tracking(run)['SB.om', ac(ay)])

dimnames(ext$output) <- list(year=ay)
dimnames(ext$change) <- list(year=ay)

dat <- data.table(model.frame(ext, drop=TRUE))


inspect  <- function(run) {

  # GET args
  spread(args(run))

  # SET years
  ays <- as.numeric(vy)
  dys <- ays - data_lag
  mys <- ays + management_lag

  # EXTRACT FLQs
  mets <- FLQuants(
    # ay, from dy
    basis=tracking(run)['metric.hcr', ac(ays)],
    # ay
    rule=tracking(run)['rule.hcr', ac(ays)],
    # ay
    decision=tracking(run)['decision.hcr', ac(ays+1)],
    # mys
    output=tracking(run)['hcr', ac(ays+1)],
    # msy
    management=tracking(run)['fwd', ac(ays+1)],
    # ay
    status=tracking(run)['SB.om', ac(mys)])
  
  dat <- data.table(model.frame(ext, drop=TRUE))

  return(dat)
}

dat <- inspect(run)



library(latex2exp)

# PLOT basis (metric) distribution by year, color by rule

# plot_decision_basis

ggplot(dat, aes(x=ISOdate(year, 1, 1), y=basis, group=year)) +
  geom_hline(yintercept=unlist(args(control(run)$hcr)[c('lim', 'bufflow', 'buffupp')]),
    linetype=2, alpha=0.5) +
  geom_violin(fill='gray', alpha=0.3) +
  geom_point(aes(colour=factor(rule)), alpha=0.4) +
  ylab("Depletion") + xlab("") + ylim(0, 1) +
  scale_color_manual(name="Decision", values=flpalette_colours(4),
    labels=setNames(TeX(c('$D > L$', '$L \\leq D < B_{low}$',
      '$B_{low} \\leq D < B_{upper}$', '$D \\geq B_{upper}$')), nm=1:4))






# buffer
ggplot(dat, aes(x=ISOdate(year, 1, 1), y=basis, group=year)) +
  geom_hline(yintercept=unlist(args(control(run)$hcr)[c('lim', 'bufflow', 'buffupp')]),
    linetype=2, alpha=0.5) +
  geom_violin(fill='gray', alpha=0.3) +
  geom_point(aes(colour=factor(rule)), alpha=0.4) +
  ylab("Depletion") + xlab("") + ylim(0, 1) +
  scale_color_manual(name="Decision", values=flpalette_colours(4),
    labels=setNames(TeX(c('$D > L$', '$L \\leq D < B_{low}$',
      '$B_{low} \\leq D < B_{upper}$', '$D \\geq B_{upper}$')), nm=1:4))

# hockeystick
ggplot(dat, aes(x=ISOdate(year, 1, 1), y=basis, group=year)) +
  geom_hline(yintercept=unlist(args(ctrl$hcr)[c('lim', 'trigger')]),
    linetype=2, alpha=0.5) +
  geom_violin(fill='gray', alpha=0.3) +
  geom_point(aes(colour=factor(rule)), alpha=0.4) +
  ylab("Depletion") + xlab("") + ylim(0, 1) +
  scale_color_manual(name="Decision", values=flpalette_colours(4),
    labels=setNames(TeX(c('$D > L$', '$L \\leq D < T$',
      '$D \\geq T$')), nm=1:3))

# PLOT rule (proportion) by year

ggplot(dat[, .N, by=.(year, rule)], aes(x=ISOdate(year,1,1), y=N / 100)) +
  geom_bar(aes(fill=factor(rule)), sta='identity') + 
  scale_fill_manual(name="Decision", values=flpalette_colours(4),
    labels=setNames(TeX(c('$D < L$', '$L \\leq D < B_{low}$',
      '$B_{low} \\leq D < B_{upper}$', '$D \\geq B_{upper}$')), nm=seq(1, 4))) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion by rule category") + xlab("")


dau <- melt(dat, id.vars=c('year', 'iter'),
  measure.vars=c("decision", "output", "management"), variable.name = "step", value.name = "data")

# decision
# output
# management

ggplot(dau, aes(x=year, y=data, group=year)) +
  geom_violin() +
  facet_grid(step~., scales='free')








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

