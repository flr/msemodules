# om.R - DESC
# msemodules/R/om.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# asFLom {{{

setGeneric("asFLom", function(object, ...) standardGeneric("asFLom"))

# FLa4a

.asFLoma4aFitSA <- function() {

setMethod("asFLom", signature(object="a4aFitSA"),

function(object, stk, sr_model, fy, it=1, nsq=3, sa2sr=TRUE,
  refpts=list(bsafe = "spr.40", blim = "spr.20", fmsy = "f0.1"), ...){

  # CHECK for FLa4a
  if (!requireNamespace("FLa4a", quietly = TRUE)) {
    stop("Package 'FLa4a' is required for this method. Please install it.",
      call. = FALSE)
  }

  # DIMS
  dy <- dims(stk)$maxyear + 1
  
  # update stock with simulation from covar
  sim <- simulate(object, it)
  stk <- stk + sim

  # S/R
  if(isTRUE(sa2sr)){
    om.sr <- fmle(as.FLSR(stk, model=sr_model), control = list(trace = 0))
  } else {
    om.sr <- fmle(as.FLSR(qapply(stk, iterMedians), model=sr_model),
      control = list(trace = 0))
  }

  # SR deviances
  om.srdevs <- append(exp(residuals(om.sr)),
    rlnormar1(it, sdlog=sqrt(yearVars((residuals(om.sr)))), years=seq(dy, fy)))

  # BRP
  om.brp <- FLBRP(stk, sr=om.sr)
  rp <- refpts(om.brp)

  # ADD refpts
  rpnms <- unique(c(dimnames(rp)$refpt, unlist(refpts)))

  if(length(rpnms) > 7){
    rp <- rp[c(1:nrow(rp), rep(1, length(rpnms) - 7))]
    dimnames(rp)$refpt <- rpnms
    refpts(om.brp) <- rp
  }
  
  om.brp <- brp(om.brp)

  # NOTE: Blim, Btrigger
  rp <- remap(refpts(om.brp),
    map=list(BLIM=c(refpts$blim, "ssb"), BSAFE = c(refpts$bsafe, "ssb"),
        FMSY = c(refpts$fmsy, "harvest"), SBMSY = c("msy", "ssb"),
        BMSY = c("msy", "biomass"), B0 = c("virgin", "biomass"),
        SB0 = c("virgin", "ssb")))

  # OM
  om <- FLom(stock=stk, refpts=rp, model=sr_model, params=params(om.sr),
    deviances=om.srdevs, name=name(stk))

  # EXTEND
  om <- fwdWindow(om, end=fy, nsq=nsq, ...)
  
  return(om)
  }
)
}

# }}}

# runOM {{{

#' Run a life-history operating model over a historical trajectory.
#'
#' Internal helper that builds an equilibrium population from `lhs`, scales its
#' fishing mortality to the supplied history, and projects it forward with
#' recruitment deviances. Fishing mortality is capped at 90% of the crash
#' reference point before projection. The returned `FLStockR` contains the
#' projected `FLStock`, its stock-recruitment model, and MSY reference points.
#'
#' @param lhs A life-history specification accepted by [lhEql()].
#' @param history A matrix or compatible object used to scale equilibrium
#'   fishing mortality by time step.
#' @param deviances Recruitment deviances. Its iteration dimension determines
#'   the number of operating-model iterations to project.
#' @param ... Additional arguments forwarded to [lhEql()]. Supplying `range`
#'   prevents the helper from deriving the `minfbar` and `maxfbar` range from
#'   the equilibrium population.
#'
#' @return An `FLStockR` object containing the historical projection, a
#'   `predictModel` stock-recruitment relationship, and `FMSY`, `SBMSY`, and
#'   `MSY` reference points.
#'
#' @keywords internal
#' @noRd
runOM <- function(lhs, history, deviances, ...) {

  # VARIABLES
  its <- dims(deviances)$iter

  # CREATE equilibrium population
  eq <- lhEql(lhs, ...)

  # ADJUST fbar to history
  fbar(eq) <- history %*% fbar(eq)[catch(eq) == max(catch(eq))]

  # BUT limited by 90% Fcrash
  fbar(eq) <- qmin(fbar(eq), refpts(eq)["crash", "harvest"] * 0.9)

  # COERCE into om FLStock
  om <- as(eq, "FLStock")

  # SET catch >= 0
  stock.n(om)[stock.n(om) < 0] <- 0.0001
  catch.n(om)[catch.n(om) < 0] <- 0.0001
  landings.n(om)[landings.n(om) < 0] <- 0.0001

  # SET fbar range
  if(!"range" %in% names(list(...)))
    range(om, c("minfbar", "maxfbar")) <- 
      mean(ages(m(om))[harvest(om) == expand(fapex(om), age=dimnames(om)$age)])

  # PROJECT for history
  om <- ffwd(propagate(om, its), sr=eq, fbar=fbar(eq)[, -1],
    deviances=deviances)

  # EXTRACT refpts
  rps <- FLPar(FMSY=fmsy(eq), SBMSY=sbmsy(eq), MSY=msy(eq))

  # OUTPUT FLStockR
  res <- FLStockR(om, sr=as(eq, "predictModel"), refpts=rps)

  return(res)
}
# }}}

# initiate {{{

#' Initialize an unfished population for a target virgin biomass.
#'
#' Calculates the unfished equilibrium numbers-at-age required for each
#' iteration to attain `B0`. The calculation applies natural mortality between
#' ages and treats the oldest age as a plus group. It then installs a
#' Beverton-Holt stock-recruitment model parameterized by steepness, virgin
#' recruitment, and spawning-stock biomass.
#'
#' The input `FLBiol` must already contain valid weight-at-age (`wt`), natural
#' mortality (`m`), spawning timing (`spwn`), and maturity-at-age (`mat`)
#' values in its first year. Only first-year abundances and stock-recruitment
#' parameters are changed.
#'
#' @param biol An `FLBiol` object whose first-year biology defines the
#'   unfished equilibrium population.
#' @param B0 Target virgin total biomass. A scalar is recycled across
#'   iterations; otherwise values are matched to the iteration dimension.
#' @param h Beverton-Holt steepness. The default is 0.75.
#'
#' @return `biol`, with first-year `n` set to the equilibrium abundance vector
#'   and `sr` parameterized with the corresponding virgin recruitment (`R0`)
#'   and spawning biomass (`v`).
#' @export
#'
#' @examples
#' data(ple4.biol)
#' \dontrun{
#' initiate(ple4.biol, B0=450000)
#' # Sets up parallel
#' # plan(multisession, workers=4)
#' initiate(propagate(ple4.biol, 100), B0=runif(100, 3e5, 5e5))
#' }

initiate <- function(biol, B0, h=0.75) {
  
  # SET iters
  dms <- dims(biol)
  its <- dms$iter
  na <- dms$age
  B0 <- rep(B0, length=its)

  # SETUP initial sr
  sr(biol) <- predictModel(model=bevholtss3()$model,
    params=propagate(FLPar(s=h, R0=NA, v=NA), its))

  # SOLVE R0 for B0, WT, M + F
  foo <- function(R0, n, m, wt, b0) {

    n[1] <- R0

    for(a in seq(2, length(n)))
      n[a] <- n[a - 1] * exp(-m[a - 1])
    # plusgroup
    n[a] <- n[a] / (1 - exp(-m[a]))

    return((c(b0) - sum(wt * n)))
  }
 
  # TODO: DEAL with iters only in biol, B0 and h

  # INITIAL value, assumes 1,000 rec per t SSB.
  init <- B0

  # RUN for iters
  res <- foreach(i=seq(its), .combine=c) %dofuture% {

    # EXTRACT to vectors
    m <- c(iter(m(biol)[, 1], i))
    wt <- c(iter(wt(biol)[, 1], i))
    n <- c(iter(n(biol)[, 1], i))

    res <- uniroot(foo, c(1, 1e12), n=n, m=m, wt=wt, b0=B0[i])
    print(i)
    res$root
  }

  # RECONSTRUCT initial population

  # rec
  n(biol)[1, 1] <- res
    
  # n
  for(a in seq(2, na))
    n(biol)[a, 1] <- n(biol)[a - 1, 1] * exp(-m(biol)[a - 1, 1])

  # pg
  n(biol)[na, 1] <- n(biol)[na, 1] / (1 - exp(-m(biol)[na, 1]))
  
  # ADD R0 param
  params(sr(biol))$R0 <- res
  params(sr(biol))$v <- ssb(biol)[,1]

  # RETURN FLBiol
  return(biol)
}
# }}}

# deplete {{{

#' Set a population to a specified depletion level.
#'
#' Builds an `FLBRP` object from first-year biological inputs and the supplied
#' catch selectivity, parameterizes it with the `FLBiol` Beverton-Holt
#' stock-recruitment relationship, and finds the fishing mortality producing
#' the requested fraction of virgin biomass. The resulting equilibrium
#' numbers-at-age replace the first-year abundance in `biol`; calculated
#' reference points are retained as its `"refpts"` attribute.
#'
#' `dep` is interpreted as a fraction of virgin total biomass, so `0` denotes
#' the lowest biomass reachable by the reference-point calculation and `1`
#' denotes the unfished biomass. For depletion calculations, multiple
#' fisheries should first be represented by a combined catch selectivity.
#'
#' @param biol An initialized `FLBiol` whose first year supplies biological
#'   inputs and a Beverton-Holt stock-recruitment relationship.
#' @param sel Catch selectivity at age. Only its first year is used.
#' @param dep Target depletion as a proportion of virgin biomass, between zero
#'   and one. Values are applied by iteration.
#' @param minfbar Youngest age included in the fishing-mortality average.
#'   Defaults to the minimum age in `biol`.
#' @param maxfbar Oldest age included in the fishing-mortality average.
#'   Defaults to the maximum age in `biol` and is capped at that age.
#'
#' @return `biol`, with first-year `n` at the equilibrium depletion level and
#'   an attached `"refpts"` attribute containing the `FLBRP` reference points.
#' @export
#'
#' @examples
#' data(ple4.biol)
#' ini <- initiate(propagate(ple4.biol, 100), B0=runif(100, 3e5, 5e5))
#' dep <- deplete(ini, sel=catch.sel(ple4)[, 10], dep=runif(100, 0.10, 0.90))

deplete <- function(biol, sel, dep, minfbar=dims(biol)$min, 
  maxfbar=dims(biol)$max) {

  # GET dims
  dm <- dim(n(biol))
  myr <- dims(biol)$minyear
  mag <- dims(biol)$max

  # CHECK inputs
  # dep ~ [0,1]
  if(!all(dep <= 1 & dep >= 0))
    stop("depletion (dep) must fall between 0 and 1.")
  # maxfbar <= mag
  if(maxfbar > mag)
    maxfbar <- mag

  # EXTRACT slots
  waa <- wt(biol)[, 1]
  maa <- m(biol)[, 1]
  mat <- mat(biol)[, 1]
  msp <- spwn(biol)[, 1]
  sr <- sr(biol)

  # USE only first year of sel
  sel <- sel[, 1]

  # FLBRP
  brp <- FLBRP(stock.wt=waa, landings.wt=waa, discards.wt=waa, bycatch.wt=waa,
    mat=mat, landings.sel=sel, m=maa,
    discards.sel=sel %=% 0, bycatch.harvest=sel %=% 0,
    harvest.spwn=maa %=% 0, m.spwn=maa %=% 0,
    availability=maa %=% 1,
    range=c('minfbar'=minfbar, 'maxfbar'=maxfbar, 'plusgroup'=mag))

  # ADD sr as bevholt(a,b)
  psr <- params(sr)
  npsr <- abPars("bevholt", spr0=psr$v / psr$R0, s=psr$s, v=psr$v)
  model(brp) <- bevholt()$model
  params(brp) <- FLPar(a=npsr$a, b=npsr$b)
  
  # ADD Btgt
  brp <- brp(brp)
  refpts(brp, "target", "biomass") <- c(refpts(brp)['virgin', 'biomass',]) * dep
  
  # SET fbar as ftarget by iter
  fbar(brp) <- FLQuant(c(refpts(brp)['target', 'harvest',]),
    dim=c(1, 1, 1, 1, 1, dm[6]))

  # COMPUTE & ASSIGN stock.n
  n(biol)[,1] <- stock.n(brp)

  # TEST: abs(range(c(tb(biol)[,1]) / (B0 * dep)))

  attr(biol, "refpts") <- refpts(brp)

  return(biol)
}
# }}}

# simulator {{{

#' Simulate an age-structured operating model.
#'
#' `simulator()` creates an initial unfished population for every iteration,
#' depletes it to the requested initial biomass fraction, and then projects it
#' with [fwd()] using the supplied fisheries and historical control. It derives
#' a combined catch selectivity when more than one fishery is present, scales
#' fleet effort to the fishing mortality associated with the depletion target,
#' and stores reference points, priors, and recruitment deviances alongside
#' the projection.
#'
#' Large iteration sets are split into blocks of 250 so that iterations can be
#' evaluated through the active `doFuture` backend. A progressr progressor is
#' advanced once per completed block. If `invalk` is supplied, the function
#' also draws 250 length samples per fishery from simulated catch-at-age.
#'
#' @param biol An `FLBiol` that supplies the age-structured biology, including
#'   weight, mortality, maturity, and spawning timing.
#' @param fisheries An `FLFisheries` object or list of fisheries defining catch
#'   selectivity, catches, and effort for the forward projection.
#' @param history Historical fishing controls. Objects other than `fwdControl`
#'   are coerced to `fwdControl` before projection.
#' @param B0 Virgin total biomass by iteration. Its length determines the
#'   number of simulated iterations.
#' @param h Beverton-Holt steepness by iteration.
#' @param dep Initial depletion as a proportion of virgin biomass by iteration.
#'   The default (`0`) starts each population at the lowest biomass calculated
#'   by [deplete()].
#' @param iter Nominal number of iterations. The current implementation derives
#'   the number of simulations from the length of `B0`; this argument is
#'   retained for compatibility.
#' @param sigmaR Standard deviation of lognormal recruitment deviances. It is
#'   also used to perturb initial abundances below the plus group.
#' @param rho Autocorrelation of generated recruitment deviances. Used only
#'   when `deviances` is not supplied.
#' @param B0change Optional annual multiplier series for virgin biomass and
#'   virgin recruitment, allowing carrying capacity to vary through time.
#' @param invalk Optional inverse age-length key passed to [lenSamples()] to
#'   generate length samples from catch-at-age.
#' @param deviances Recruitment deviances for the projection. By default these
#'   are generated with [ar1rlnorm()] over the biological object's years.
#' @param minfbar Youngest age included in the fishing-mortality average used
#'   while calculating depletion.
#' @param maxfbar Oldest age included in the fishing-mortality average used
#'   while calculating depletion.
#'
#' @return A list containing the simulated `biol`, `fisheries`, iteration-level
#'   `priors` (`B0`, `dep`, and `h`), recruitment `deviances`, and `refpts`.
#'   When `invalk` is supplied, it also contains simulated catch length samples
#'   in `lengths`.
#' @export
#' @examples
#' NULL

simulator <- function(biol, fisheries, history, B0, h, dep=0,
  iter=dims(biol)$iter, sigmaR=0, rho=0, B0change=NULL, invalk=NULL,
  deviances=ar1rlnorm(rho=rho, years=dimnames(biol)$year,
    meanlog=0, sdlog=sigmaR),
  minfbar=dims(biol)$min, maxfbar=dims(biol)$max) {

  # DIMS
  nage <- dims(biol)$age
  nyr <- dims(biol)$year
  nfs <- length(fisheries)
  its <- length(c(B0))
  bls <- split(seq(its), ceiling(seq_along(seq(its)) / 250))

  # CHECK match history and first year with N
  
  # SET progresssor
  p <- progressor(length(bls))

  # PROPAGATE objects, if needed

  if(dims(biol)$iter == 1)
    biols <- lapply(bls, function(x) propagate(biol, length(x)))
  else if(dims(biol)$iter == its)
    biols <- lapply(bls, iter, obj=biol)
  else
    biols <- list(`1`=biol)

  if(dims(fisheries[[1]])$iter == 1)
    fisheries <- lapply(bls, function(x) lapply(fisheries,
      function(y) propagate(y, iter=length(x))))
  else if(dims(fisheries[[1]])$iter == its)
    fisheries <- lapply(bls, iter, obj=fisheries)
  else
    fisheries <- list(`1`=fisheries)

  # LOOP over iter blocks

  sim <- foreach(i=names(bls), .combine=.lcombine,
    .multicombine=TRUE) %dofuture% {

    # SET iters
    it <- bls[[i]]
    ni <- length(it)

    # GET objects by iter group from lists
    bio <- biols[[i]]
    fis <- fisheries[[i]]
    
    # INITIATE N0
    nbio <- initiate(bio, B0=B0[it], h=h[it])

    # COMBINED selectivity for depletion, catch 1 across fisheries
    if (nfs > 1)
    sel <- Reduce("+", lapply(fis, function(x) catch.sel(x[[1]])) *
      lapply(fis, function(x) catch.n(x[[1]]))) /
      Reduce("+", lapply(fis, function(x) catch.n(x[[1]])))
    else
      sel <- catch.sel(fis[[1]][[1]])

    # RESCALE selex to 1
    sel <- sel %/% apply(sel, 2:6, max)

    # DEPLETE to dep level by iter
    nbio <- deplete(nbio, sel=sel, dep=dep[it], minfbar=minfbar, 
      maxfbar=maxfbar)

    # ADD initial age devs, no bias correction needed
    n(nbio)[-nage, 1] <- n(nbio)[-nage, 1] * rlnorm(nage - 1, 0, sigmaR)

    # MAP refpts & keep target F
    rps <- remap(nbio@refpts, MSY=c("msy", "yield"))
    targetf <- c(nbio@refpts["target", "harvest"])

    # ALTER SRR if B0(K) changes in time
  
    if(!is.null(B0change)) {

      # EXPAND FLPar
      pas <- params(sr(nbio))
      pay <- FLPar(NA, dimnames=list(params=c("s", "R0", "v"),
       year=dimnames(nbio)$year, iter=it))

      # ASSIGN first year,
      pay[, 1,]<- pas
      # steepness
      pay['s', ,] <- pas['s',]
      # SET B0 and R0 trends
      pay['v', ,] <- apply(pas$v, 2, '*', c(B0change))
      pay['R0', ,] <- apply(pas$R0, 2, '*', c(B0change))

      params(sr(nbio)) <- pay
  
      # RESCALE refpts by year
      rpy <- FLPar(NA, dimnames=list(param=dimnames(rps)$param, 
        year=dimnames(nbio)$year, iter=it))
      rpy[,1,] <- rps

      # F is constant
      rpy[c('FMSY'),]  <- rps[c('FMSY'),] 

      # BIOMASS is scaled
      for(i in c('SBMSY', 'BMSY', 'B0', 'SB0', 'MSY'))
        rpy[i,] <- apply(rps[i, ], 2, '*', c(B0change))
    }

    # SET effort to match F target
    for(i in seq(fis)) {
      effort(fis[[i]])[] <- abs(targetf)
    }

    # CONVERT history
    # TODO: DEAL w/ iters in history
    if(!is(history, "fwdControl")) {
      history <- as(history, "fwdControl")
    }
    
    # FWD w/history
    res <- suppressWarnings(fwd(nbio, fis, control=iter(history, it),
      deviances=iter(deviances, it), effort_max=1e6))

    # LEN samples
    if(!is.null(invalk)) {
      cafs <- lapply(res$fisheries, function(x) catch.n(x[[1]]) + 1e-6)
      res$lengths <- lapply(cafs, lenSamples, invALK=invalk, n=250)
    }

    # OUTPUT
    res$priors <- data.table(B0=B0, dep=dep, h=h)
    deviances(res$biol) <- iter(deviances, it)
    res$deviances <- iter(deviances, it)
    res$refpts <- rps

    # UPDATE progressr report
    p()

    return(res)
  }

  return(sim)
}

# }}}

# .bcombine including refpts attr {{{
#' Combine biological objects while retaining reference points.
#'
#' This internal reducer combines its inputs with [combine()] and combines
#' their `"refpts"` slots separately, since the standard combination does not
#' preserve that attribute.
#'
#' @param x,y Biological simulation objects to combine.
#' @param ... Additional biological simulation objects to combine.
#'
#' @return The combined object, with a combined `"refpts"` attribute.
#'
#' @keywords internal
#' @noRd
.bcombine <- function(x, y, ...) {

  args <- c(list(x, y), list(...))
  res <- do.call(combine, args)

  attr(res, "refpts") <- Reduce(combine, lapply(args, slot, "refpts"))

  return(res)
}
# }}}

# .lcombine simulator list {{{
#' Combine blocks returned by [simulator()].
#'
#' This internal `foreach` reducer merges the biological population, fisheries,
#' priors, recruitment deviances, and reference points returned for independent
#' iteration blocks. Length samples are merged when present.
#'
#' @param x,y Simulation-result lists to combine.
#' @param ... Additional simulation-result lists to combine.
#'
#' @return A simulation-result list with all iteration blocks combined.
#'
#' @keywords internal
#' @noRd
.lcombine <- function(x, y, ...) {

  args <- c(list(x, y), list(...))

  out <- list(biol=do.call(.bcombine, lapply(args, '[[', 'biol')),
    fisheries=do.call(combine, lapply(args, '[[', 'fisheries')),
    priors=do.call(rbind, lapply(args, '[[', 'priors')),
    deviances=do.call(combine, lapply(args, '[[', 'deviances')),
    refpts=do.call(combine, lapply(args, '[[', 'refpts'))
  )

  if("lengths" %in% names(x))
    out$lengths <- do.call(combine, lapply(args, '[[', 'lengths'))

  return(out)
}
# }}}
