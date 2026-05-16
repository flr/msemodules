# isys.R - DESC
# /home/mosqu003/Active/mse_FLR/msemodules/R/isys.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


globalVariables(c("metric"))

# bank_borrow.is {{{

#' Banking and borrowing implementation system method
#'
#' Adjusts a TAC recommendation from an HCR by banking a proportion of any
#' increase or borrowing against a future quota when a decrease is large,
#' smoothing inter-annual TAC variability while respecting stock health
#' conditions.
#'
#' @param stk An [FLCore::FLStock-class] object representing the stock.
#' @param ctrl An [FLasher::fwdControl-class] object carrying the TAC
#'   recommendation to be adjusted.
#' @param args A list of management-cycle dimensionality arguments supplied
#'   automatically by [mse::mp()]. Must contain at least `iy`, `ay`, `frq`,
#'   and `management_lag`.
#' @param split Optional split-control input. Currently unused; reserved for
#'   future support of multiple [mse::mseCtrl-class] objects within a single
#'   `isys` step. Defaults to `NULL`.
#' @param rate Numeric in \eqn{[0, 1]}. Proportion of the TAC to bank or
#'   borrow when the threshold condition is met. Must be supplied; there is no
#'   default.
#' @param diff Numeric. Minimum relative change in TAC (compared to the
#'   previous HCR decision `pre`) required to trigger a banking or borrowing
#'   action. Defaults to `0.15`, i.e. a 15% change.
#' @param healthy Numeric. Minimum value of the `rule.hcr` tracking metric
#'   required for the stock to be considered healthy enough for banking or
#'   borrowing to be applied. Defaults to `1`.
#' @param tracking An [FLCore::FLQuant-class] used to record intermediate
#'   values and decisions during MP evaluation, passed through and updated by
#'   [mse::mp()].
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`ctrl`}{The input [FLasher::fwdControl-class] object with the
#'     adjusted TAC stored in `ctrl$value`.}
#'   \item{`tracking`}{The updated tracking object, with `borrowing.isys` and
#'     `banking.isys` metrics written for the assessment year `ay`.}
#' }
#'
#' @details
#' In the initial year (`ay == iy`) the tracking metrics `borrowing.isys` and
#' `banking.isys` are initialized to zero and the reference TAC (`pre`) is set
#' to the realized catch in the year preceding the management lag.
#'
#' In subsequent years `pre` is taken from the `hcr` tracking metric for the
#' current year, and the TAC carried in `ctrl` is first corrected for any
#' amount borrowed or banked in the previous cycle before new banking/borrowing
#' is evaluated.
#'
#' Banking and borrowing are triggered as follows:
#' \itemize{
#'   \item **Borrowing**: if the (corrected) TAC falls more than `diff` below
#'     `pre` \emph{and} the stock is healthy, `rate * tac` is added to the
#'     current TAC and recorded in `borrowing.isys`.
#'   \item **Banking**: if the (corrected) TAC rises more than `diff` above
#'     `pre` \emph{and} the stock is healthy, `rate * tac` is subtracted from
#'     the current TAC and recorded in `banking.isys`.
#' }
#'
#' Only annual management cycles (`frq == 1`) are currently supported; a
#' non-annual `frq` raises an error.
#'
#' @seealso [mse::mp()], [mse::mpCtrl()], [mse::mseCtrl()]
#'
#' @examples
#' \dontrun{
#' data(sol274)
#'
#' # mpCtrl combining a hockey-stick HCR with 10% banking/borrowing
#' ctrl <- mpCtrl(
#'   est   = mseCtrl(method = perfect.sa),
#'   hcr   = mseCtrl(method = hockeystick.hcr,
#'             args = list(metric = "ssb", trigger = 42000,
#'                         output = "catch", target = 11000)),
#'   isys  = mseCtrl(method = bank_borrow.is,
#'             args = list(rate = 0.10, healthy = 2, diff = 0.05)))
#'
#' run <- mp(om, control = ctrl, args = list(iy = 2021, fy = 2035))
#'
#' plot(om, list(BaB = run))
#'
#' # Inspect the TAC-setting steps from the tracking object
#' items <- c("year", "hcr", "banking.isys", "borrowing.isys", "isys", "fwd")
#' dcast(
#'   tracking(run)[metric %in% items[-1],
#'     .(data = mean(data)), by = .(year, metric)],
#'   year ~ metric, value.var = "data")[, ..items]
#' }
#'
#' @author Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#' @keywords manip

bank_borrow.is <- function(stk, ctrl, args, split=NULL, rate = NULL, diff = 0.15,
  healthy=1, tracking) {

  # DIMS
  spread(args)

  # STOP if frq > 1
  if(frq > 1)
    stop("Banking & borrowing currently assumes annual management (frq=1)")

  # GET TAC
  tac <- ctrl$value

  # CORRECT for previous borrowing or banking
  if (ay == iy) {

    # START tracking
    track(tracking, "borrowing.isys", year=ay) <- 0
    track(tracking, "banking.isys", year=ay) <- 0

    # SET initial
    pre <- c(areaSums(unitSums(seasonSums(window(catch(stk),
      start = ay - management_lag, end = ay - management_lag)))))
  
  } else {

    # GET this year's HCR decision
    pre <- tracking[metric == 'hcr' & year == ay, data]

    # GET banked or borrowed amounts
    borrowed <- tracking[metric == 'borrowing.isys' & year == ay - frq, data]
    banked <- tracking[metric == 'banking.isys' & year == ay - frq, data]

    # CORREDCT tac for ongoing banking & borrowing
    tac <- tac - borrowed + banked
  }

  # CHECK status
  id <- tracking[metric == "rule.hcr" & year == ay, data] >= healthy

  # IF lower TAC, THEN borrow
  borrowing <- ifelse(tac < pre * (1 - diff) & id, tac * rate, 0)

  # TRACK amount being borrowed from mys
  track(tracking, "borrowing.isys", ay) <- borrowing

  # IF higher TAC, THEN bank
  banking <- ifelse(tac > pre * (1 + diff) & id, tac * rate, 0)
  
  # TRACK amount being banked into mys
  track(tracking, "banking.isys", ay) <- banking

  # CORRECT TAC
  ctrl$value <- tac + borrowing - banking

  # split TODO: isys to accept multiple mseCtrls
  # res <- split.is(stk, ctrl, split, quant = "catch", args, tracking)
  # return(list(ctrl = res$ctrl, tracking = res$tracking))

  return(list(ctrl = ctrl, tracking = tracking))
}
# }}}

# effort.is {{{

#' Effort-multiplier implementation system method
#'
#' Converts a fishing-mortality target set by an HCR into a relative effort
#' multiplier by scaling against a reference fishing mortality computed from
#' recent years. Passing a relative effort target to [FLasher::fwd()] rather
#' than an absolute F smooths the response to inter-annual variability in
#' fleet dynamics.
#'
#' @param stk An [FLCore::FLStock-class] object containing stock and fishery
#'   information.
#' @param ctrl An [FLasher::fwdControl-class] object whose `quant` slot must be
#'   one of `"f"`, `"fbar"`, or `"effort"` and whose `value` slot carries the
#'   absolute F target from the HCR.
#' @param nyears Integer. Number of years over which the mean reference fishing
#'   mortality is calculated. Defaults to `args$nsqy` as supplied by
#'   [mse::mp()].
#' @param Fref Numeric or [FLCore::FLQuant-class]. Reference fishing mortality
#'   used as the denominator when computing the multiplier. Defaults to the
#'   mean of [FLCore::fbar()] over the `nyears` years ending at `dy` (the
#'   final data year).
#' @param args A list of management-cycle dimensionality arguments supplied
#'   automatically by [mse::mp()], including at least `dy` and `nsqy`.
#' @param tracking An [FLCore::FLQuant-class] used to record intermediate
#'   values and decisions during MP evaluation, passed through unchanged by
#'   this method.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`ctrl`}{The input [FLasher::fwdControl-class] object with `value`
#'     replaced by the dimensionless effort multiplier and `relYear` set to
#'     `ctrl$year - data_lag` so that [FLasher::fwd()] interprets the target
#'     relative to the correct historical year.}
#'   \item{`tracking`}{The tracking object, returned unchanged.}
#' }
#'
#' @details
#' The effort multiplier is defined as
#' \deqn{m = \frac{Ftarget}{F_ref}}{mult = target / Fref}
#' where \eqn{F_target} is the value in `ctrl$value` as set by the
#' HCR and \eqn{F_ref} is the mean [FLCore::fbar()] over the most
#' recent `nyears` years of data.
#'
#' Setting `ctrl$relYear` tells [FLasher::fwd()] to express the effort target
#' relative to the effort observed in that reference year, effectively turning
#' the absolute F recommendation into a proportional change instruction.
#'
#' An error is raised immediately if `ctrl$quant` is not one of `"f"`,
#' `"fbar"`, or `"effort"`.
#'
#' @seealso [mse::mp()], [mse::split.is()], [mse::mpCtrl()], [mse::mseCtrl()],
#'   [FLCore::fbar()]
#'
#' @examples
#' \dontrun{
#' data(plesim)
#'
#' # mpCtrl combining a hockey-stick HCR (F target) with effort.is
#' ctrl <- mpCtrl(
#'   est  = mseCtrl(method = perfect.sa),
#'   hcr  = mseCtrl(method = hockeystick.hcr,
#'            args = list(metric = "ssb", trigger = 45000,
#'                        output = "fbar", target = 0.27)),
#'   isys = mseCtrl(method = effort.is, args = list(nyears = 3)))
#'
#' run      <- mp(om, control = ctrl,      args = list(iy = 2021, fy = 2035))
#' run_nois <- mp(om, control = ctrl[-3],  args = list(iy = 2021, fy = 2035))
#'
#' # Compare with and without the effort.is buffer effect
#' plot(om, effort.is = run, no_is = run_nois)
#' }
#'
#' @author Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#' @keywords manip

effort.is <- function(stk, ctrl, nyears=args$nsqy, 
  Fref=yearMeans(fbar(stk)[, ac(seq(dy - nyears + 1, dy))]), 
  args, tracking) {

  # CHECK ctrl sets F or effort
  if(!ctrl$quant %in% c("fbar", "f", "effort"))
    stop("'effort.is' can only accept ctrl set on 'f', 'fbar' or 'effort'")

  # EXTRACT args
  spread(args)

	# target to reach defined by HCR
	trgt <- ctrl$value
  
  # multiplier
	mult <- trgt / c(Fref)

  # new control file, in relative terms
  ctrl$value <- mult
  ctrl$relYear <- ctrl$year - data_lag

  # TODO: TEST and SET for relative ctrl$fishery, biol, catch

  return(list(ctrl = ctrl, tracking = tracking))
}
# }}}
