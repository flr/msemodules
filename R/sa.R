# sa.R - Status and trend estimators and indicators.
# msemodules/R/sa.R

# Copyright European Union, 2018-2021
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

globalVariables(c("y0", "dy"))

# shortcut.sa {{{

#' Shortcut stock assessment including error on a status metric
#'
#' Performs a simplified or "shortcut" stock assessment by truncating an
#' operating model stock object to the current data year, computing a selected
#' stock metric, applying stochastic deviations, and returning the resulting
#' indicator together with the truncated stock and updated tracking object.
#'
#' This function is intended for use inside management procedure workflows where
#' a full assessment model is not run, but an index or assessment quantity is
#' approximated from the operating model using a direct transformation of stock
#' state and observation-error-like deviations.
#'
#' @param stk An `FLStock`-like object representing the stock available to the
#'   shortcut assessment. The object is truncated to the terminal data year
#'   `dy`.
#' @param idx An index object. Currently included for interface compatibility
#'   with other assessment functions, but not used directly inside this
#'   function.
#' @param metric A character string naming the function used to compute the
#'   assessment metric from `stk`. Defaults to `"ssb"`. The named function is
#'   called via [base::do.call()].
#' @param SSBdevs An object containing deviations to apply to the computed
#'   metric. Defaults to `met %=% 1`, i.e. a unit deviation structure matching
#'   the dimensions of the computed metric. Typically this is an `FLQuant` or
#'   compatible object.
#' @param devs Deprecated alias or shorthand for `SSBdevs`. Defaults to
#'   `SSBdevs`.
#' @param args A list-like object containing dimensions and assessment settings
#'   needed by [spread()]. In particular, this function assumes that variables
#'   such as `dy`, `y0`, `ay`, and `stock` become available after calling
#'   `spread(args)`.
#' @param tracking A tracking object updated in-place through [track()] to record
#'   convergence of the shortcut assessment.
#' @param ... Additional arguments passed to the metric function named in
#'   `metric`.
#'
#' @return
#' A named list with components:
#' \itemize{
#'   \item `stk`: the stock object truncated to the terminal data year `dy`;
#'   \item `ind`: an `FLQuants` object containing the derived indicator series,
#'     named according to `metric`;
#'   \item `tracking`: the updated tracking object with convergence information
#'     recorded.
#' }
#'
#' @details
#' The function proceeds as follows:
#' \enumerate{
#'   \item values in `args` are unpacked using [spread()];
#'   \item the stock object is truncated to `end=dy` using [window()];
#'   \item the requested metric is computed by calling the function named in
#'     `metric` on `stk`;
#'   \item the result is collapsed across units using [unitSums()];
#'   \item stochastic deviations are applied over the period from `y0` to `dy`;
#'   \item the resulting quantity is wrapped in an [FLQuants()] object and named
#'     with the selected metric;
#'   \item convergence is recorded in `tracking` by setting
#'     `"conv.est"` to `1` for year `ay` and biological unit `stock`.
#' }
#'
#' @section Assumptions:
#' \itemize{
#'   \item the function named in `metric` exists and accepts `stk` as its first
#'     argument;
#'   \item `devs` or `SSBdevs` are dimensionally compatible with the metric
#'     computed from `stk`;
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item `idx` is currently unused but may be retained for consistency with
#'     other assessment-method signatures.
#'   \item `devs` is currently used to build the indicator, even though
#'     `SSBdevs` is the more explicit argument name.
#' }
#'
#' @seealso [shortcut_devs()], [FLCore::ssb()]
#'
#' @examples
#' \dontrun{
#' # dataset contains both OM (FLom) and OEM (FLoem)
#' data(plesim)
#' # Create shortcut deviances for F and SSB
#' devs <- shortcut_devs(om, Fcv = 0.3, Fphi = 0.6)
#' # Set control: sa and hcr
#' control <- mpCtrl(list(
#'   est = mseCtrl(method=shortcut.sa,
#'    args=list(devs=devs$SSB)),
#'   hcr = mseCtrl(method=hockeystick.hcr, args=list(lim=0,
#'   trigger=14000, target=0.18))))
#' # Runs mp
#' tes <- mp(om, ctrl=control, args=list(iy=2021, fy=2026))
#' }
#'
#' @author Iago Mosqueira, WMR; Ernesto Jardim, IPMA
#' @keywords methods

shortcut.sa <- function(stk, idx, metric="ssb", SSBdevs=met %=% 1, devs=SSBdevs,
  args, tracking, ...) {

  # DIMS
  spread(args)

  # SUBSET oem stock
  stk <- window(stk, end=dy)

  # COMPUTE 'metric'
  met <- unitSums(do.call(metric, c(list(stk), list(...))))

  # CREATE ind from metric and deviances
  ind <- FLQuants(met * window(devs, start=y0, end=dy))

  # NAME as metric
  names(ind) <- metric

  # TRACK 'convergence'
  track(tracking, "conv.est", year=ay, biol=stock) <- 1
 
  list(stk=stk, ind=ind, tracking=tracking)
}

# }}}

# shortcut_devs {{{

#' Generate deviation series for shortcut assessments
#'
#' Creates stochastic multiplicative deviation series for use in shortcut stock
#' assessments, typically representing uncertainty or observation error applied
#' to fishing mortality and spawning stock biomass metrics.
#'
#' The deviations are generated as lognormal autoregressive series and returned
#' as an `FLQuants` object with components for `F` and `SSB`.
#'
#' @param om An operating model object used only to define dimensions for the
#'   generated deviation series. The function uses `dims(om)$iter` for the
#'   number of iterations and `dimnames(om)$year` for the year dimension.
#' @param Fcv Numeric. Standard deviation on the log scale for the fishing
#'   mortality deviation process. Defaults to `0.212`.
#' @param Fphi Numeric. Autocorrelation parameter for the fishing mortality
#'   deviation process. Defaults to `0.423`.
#' @param SSBcv Numeric. Standard deviation on the log scale for the spawning
#'   stock biomass deviation process. Defaults to `0`.
#' @param SSBphi Numeric. Autocorrelation parameter for the spawning stock
#'   biomass deviation process. Currently not used in the implementation, which
#'   fixes the `SSB` autocorrelation at `0`. Defaults to `0`.
#' @param bias.correct Logical. If `TRUE`, applies bias correction when
#'   generating the lognormal AR(1) deviates. Passed to [rlnormar1()].
#'   Defaults to `FALSE`.
#'
#' @return
#' An `FLQuants` object with two elements:
#' \itemize{
#'   \item `F`: a lognormal AR(1) deviation series for fishing mortality;
#'   \item `SSB`: a lognormal deviation series for spawning stock biomass.
#' }
#'
#' @details
#' The function generates deviations using [rlnormar1()], with:
#' \itemize{
#'   \item `n = dims(om)$iter` iterations;
#'   \item `years = dimnames(om)$year`;
#'   \item `meanlog = 0`, so the process is centered on multiplicative value 1
#'     on the log scale before optional bias correction;
#'   \item separate variance and autocorrelation settings for `F` and `SSB`.
#' }
#'
#' The returned object is suitable for direct use in [shortcut.sa()], where the
#' deviations can be applied to derived stock metrics.
#'
#' @section Current implementation note:
#' Although `SSBphi` is provided as an argument, the function currently hardcodes
#' `rho=0` for the `SSB` process. This means the SSB deviations are independent
#' through time regardless of the value supplied to `SSBphi`.
#'
#' @seealso [shortcut.sa()], [rlnormar1()], [FLCore::FLQuants()]
#'
#' @examples
#' \dontrun{
#' # Generate default shortcut deviations
#' devs <- shortcut_devs(om)
#'
#' # More variable and autocorrelated F deviations
#' devs <- shortcut_devs(om, Fcv = 0.3, Fphi = 0.6)
#'
#' # With bias correction
#' devs <- shortcut_devs(om, bias.correct = TRUE)
#' }
#'
#' @author Iago Mosqueira, WMR; Ernesto Jardim, IPMA
#' @keywords models
#' @export

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
