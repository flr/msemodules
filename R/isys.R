# isys.R - DESC
# /home/mosqu003/Active/mse_FLR/msemodules/R/isys.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# bank_borrow.is {{{

#' bank_borrow.isys
#'
#' Banking and Borrowing Mechanism for TAC Adjustment
#'
#' @param stk `FLStock` object representing the stock dat.
#' @param ctrl `fwdControl` object for forward projection controls.
#' @param args List of additional arguments required for the execution.
#' @param split Optional. Input for splitting controls (default is `NULL`).
#' @param rate Numeric. Proportion rate for borrowing and banking adjustments (default is `NULL`).
#' @param diff Numeric. Tolerance difference for determining borrowing or banking actions (default is 0.15).
#' @param healthy Numeric. Threshold for determining stock health (default is 3).
#' @param tracking Data structure. Object to track borrowing and banking values over time.
#'
#' @return A list containing:
#' \item{ctrl}{Modified `fwdControl` object with updated TAC values.}
#' \item{tracking}{Updated tracking data structure with recorded borrowing and banking amounts.}
#'
#' @details
#' The function performs the following:
#' - Checks if management assumes annual (frq = 1). Raises an error if not.
#' - Retrieves and adjusts TAC for previous banking and borrowing actions.
#' - Verifies the stock's health status (`healthy`) using tracking data to determine eligibility for banking/borrowing.
#' - Calculates borrowing if TAC is significantly lower than the reference (`pre`) and banking if TAC is significantly higher.
#' - Updates the `ctrl` and `tracking` objects with the adjusted TAC and banking/borrowing amounts.
#'
#' @note
#' Banking and borrowing adjustments are applied only when the frequency (`frq`) is equal to 1.
#' @seealso
#'
#' @examples
#' # Example dataset
#' data(sol274)
#' 
#' # Sets up an mpCtrl using hockeystick(catch~ssb) + bank_borrow(10%)
#' ctrl <- mpCtrl(est = mseCtrl(method=perfect.sa),
#'   hcr = mseCtrl(method=hockeystick.hcr, args=list(metric="ssb", trigger=42000, 
#'     output="catch", target=11000)),
#'   isys = mseCtrl(method=bank_borrow.is, args=list(rate=0.10, healthy=2, diff=0.05)))
#' 
#' # Runs the MP
#' run <- mp(om, control=ctrl, args=list(iy=2021, fy=2035))
#' 
#' # Plot time series
#' plot(om, list(BaB=run))
#' 
#' # Observe from tracking the TAC-setting steps
#' items <- c("year", "hcr", "banking.isys", "borrowing.isys", "isys", "fwd")
#' 
#' dcast(tracking(run)[metric %in% items[-1], .(data=mean(data)), by=.(year, metric)],
#'   year~metric, value.var='data')[, ..items]

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
    track(tracking, "borrowing.isys", ac(ay)) <- 0
    track(tracking, "banking.isys", ac(ay)) <- 0

    # SET initial
    pre <- c(areaSums(unitSums(seasonSums(window(catch(stk),
      start = ay - management_lag, end = ay - management_lag)))))
  
  } else {

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

#' Calculate Effort Multiplier for Harvest Control Rules
#'
#' Converts a harvest control rule (HCR) fihsing mortality target into a multiplier
#' relative to a (dynamic) reference fishing mortality level.
#'
#' @param stk An FLStock object containing stock information
#' @param ctrl A fwdControl object setting an F or effort target.
#' @param Fref Numeric or FLQuant. Reference fishing mortality for scaling. Defaults to
#'   the mean of fbar over the specified years.
#' @param nyears Integer. Number of years to use for calculating mean Fref.
#'   Default is taken from `args$nsqy`
#' @param args A list containing dimensionality arguments, passed on by mp().
#' @param tracking An FLQuant used for tracking indicators, intermediate values, and decisions during MP evaluation.
#'
#' @return A list containing:
#'   \describe{
#'     \item{ctrl}{Updated control list with effort multiplier in `value` 
#'       and adjusted years in `rel.year`}
#'     \item{tracking}{Updated tracking object}
#'   }
#'
#' @details
#' The function calculates the effort multiplier as:
#' \deqn{mult = \frac{target}{F_{ref}}}{mult = target / Fref}
#'
#' This multiplier is then stored in the control object with years
#' adjusted backwards by `data_lag` to account for assessment lag.
#'
#' @note
#' The `ctrl$quant` parameter must be one of "f", "fbar", or "effort".
#' Other values will cause the function to stop with an error.
#'
#' @seealso
#'   \code{\link{fbar}} for calculating mean fishing mortality,
#' @examples
#' # Example dataset
#' data(plesim)
#' 
#' # Sets up an mpCtrl using hockeystick(fbar~ssb)
#' ctrl <- mpCtrl(
#'   est = mseCtrl(method=perfect.sa),
#'   hcr = mseCtrl(method=hockeystick.hcr, args=list(metric="ssb", trigger=45000, 
#'     output="fbar", target=0.27)),
#'   isys = mseCtrl(method=effort.is, args=list(nyears=3)))
#' 
#' # Runs mp between 2021 and 2035
#' run <- mp(om, control=ctrl, args=list(iy=2021, fy=2035))
#' 
#' # Runs mp without effort.is 'nyears' buffer effect
#' run_nois <- mp(om, control=ctrl[-3], args=list(iy=2021, fy=2035))
#' 
#' # Plots results
#' plot(om, effort.is=run, no_is=run_nois)

effort.is <- function(stk, ctrl, Fref=yearMeans(fbar(stk)[, ac(seq(dy - nyears, dy))]), 
  nyears=args$nsqy, args, tracking) {

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
  ctrl$rel.year <- ctrl$year - data_lag

  # TODO: TEST and SET for relative ctrl$fishery, biol, catch

  return(list(ctrl = ctrl, tracking = tracking))

} # }}}
