# tracking.R - DESC
# /home/mosqu003/Active/mse_FLR/msemodules/R/tracking.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


globalVariables(c("..metord"))

# medmad {{{

medmad <- function(x) paste0(format(median(x), digits=3), " (",
  format(mad(x), digits=3) , ")")

# }}}

# inspect {{{

#' Inspect MSE Tracking Data
#'
#' The `inspect` function turns the long `tracking` table into a year a metric summary
#' on which to track the inoputs and outputs of all steps inside a call to the `mp`
#' function.
#'
#' @param tab An `FLms+e` object or a `data.table` containing tracking data. If an `FLmse` object is supplied, the `tracking` slot is extracted.
#' @param metrics A character vector specifying the metrics to include in the output. If a single character value is given and it does not match excatly any of the contained metrics, it is used to subset using the datatable::%ilike% function. The special value `"decisions"` can be used to select all metrics from `"hcr"` onward. If `NULL`, the default, all metrics are returned.
#' @param summary A function (such as `mean` or `median`) to summarize the `data` column within the tracking `data.table` across rthe `iter` dimension. Defaults to `medmad` which returns a string with `"Median (Median Absolute Deviation)"`.
#' @details
#' The function processes the table on the `tracking` slot by:
#'   - Arranging metrics based on the order in which they are produced inside `mp()`.
#'   - Subsets metrics based on user input or predefined criteria like `"decisions"`.
#'   - Aggregates the `data` column using the `summary` function provided.
#'   - Reshapes the processed `data.table` into a year-by-metrics format for easy reading.
#'
#' The predefined metric order is that of the various steps inside `mp()`: `"om"`, `"obs"`, `"est"`, `"ind"`, `"phcr"`, `"hcr"`, `"isys"`, `"tm"`, `"iem"`, `"fb"`, `"fwd"`. Extra tracks produced inside any module are placed before the one named after
#' the module that produced it.
#'
#' @return A `data.table` where rows correspond to years and columns to selected metrics, with the tracking `data` aggregated using the specified summary function.
#'
#' @examples
#' # Load example MP run
#' data(mserun)
#' # Inspect all metrics
#' inspect(run)
#' # Inspect with different summary function
#' inspect(tracking(run), summary=mean)
#' # Inspect metrics with "SB" or "sb" in their names
#' inspect(tracking(run), metrics="SB")
#' # Inspect specific metrics
#' inspect(tracking(run), metrics=c("C.om", "C.obs", "C.est"))
#' # Inspect only decision-related metrics from 'hcr'
#' inspect(tracking(run), metrics="decisions")

inspect <- function(tab, metrics=NULL, summary=medmad) {

  # EXTRACT tracking data.table
  if(is(tab, "FLmse"))
    tab <- tracking(tab)

  # SET order of metrics
  allmetrics <- tab[, unique(metric)]
 
  # Order of mp() flow
  ord <- c("om", "obs", "est", "ind", "phcr", "hcr", "isys", "tm", "iem", "fb", "fwd")

  # SUBSET and REORDER existing
  mets <- lapply(ord, function(x) allmetrics[allmetrics %like% x])

  metord <- unique(c("year", unlist(lapply(mets, function(x)
    c(x[!x %in% ord], x[x %in% ord])))))

  # SUBSET by metrics
  if(!is.null(metrics)) {
    if(identical(metrics, "decisions")) {
      metord <- metord[seq(which(metord == "hcr"), length(metord))]
    # USE %like% if no exact match
    } else if(!all(metrics %in% metord)) {
      metord <- metord[metord %ilike% tolower(metrics)]
    # EXACT match
    } else {
      metord <- metord[metord %in% metrics]
    }
  }

  # RESHAPE and REORDER
  res <- dcast(tab[metric %in% metord, .(data=summary(data)), by=.(year, metric)],
    year ~ metric, value.var='data')[, ..metord]

  return(res)
}

# }}}

# trackingFLQuants {{{

#' Convert Performance Data to FLQuant Format
#'
#' Transforms performance metrics from `FLmse` or `FLmsea` objects, stored as a
#' data.table into an `FLQuant` pr `FLQuants` object.
#'
#' @param x An object of class `FLmse`, `FLmses`, or a `data.table`.
#'
#' @return
#'   If a single management procedure (mp) is present, returns a single `FLQuant`
#'   object. If multiple management procedures are present, returns a list of
#'   `FLQuant` objects, one per mp. When `biol` column is present, it is renamed
#'   to `unit` in the resulting `FLQuant` or `FLQuants`.
#'
#' @details
#'   If `x` is of class `FLmse` or `FLmses`, the performance slot is extracted
#'   automatically. The data.table must contain columns: `statistic`, `year`, `iter`,
#'   `data`, and `mp`. Optionally, a `biol` column can be present, as in the
#'   performance table obtained from callin `mp()` on an `FLombf` OM. Values for each 
#'   `biol` will be separated using the 'unit' dimension.
#'
#'   The first (`quant`) dimension of the `FLQuant` or `FLQuants` contain the
#'    statistics, named as in the `statistic` column of the performance table.
#'    
#'    When results from multiple management procedures are present, the output is an
#'   `FLQuants` list of `FLQuant` objects, one per mp.
#'
#' @examples
#' \dontrun{
#'   # Using a data.table directly
#'   perf_dt <- data.table(
#'     statistic = rep("SSB", 4),
#'     year = rep(2025:2026, 2),
#'     iter = rep(1:2, each = 2),
#'     data = runif(4, 1000, 3000),
#'     mp = "MP1"
#'   )
#'   quant <- trackingFLQuant(perf_dt)
#'   quants <- trackingFLQuants(perf_dt)
#' }
#'
#' @seealso [performance()], [FLmse-class], [FLmses-class]
#'
#' @keywords manip

trackingFLQuant <- function(x) {
  return(trackingFLQuants(x))
}

#' @rdname trackingFLQuant

trackingFLQuants <- function(x) {

  # USE poor mans' dispatch
  if(is(x, "FLmse")) {
    res <- .coercetrackingFLQuant(copy(tracking(x)))
  } else if(is(x, "FLmses")) {
    res <- FLQuants(lapply(x, function(x)
      .coercetrackingFLQuant(copy(tracking(x)))))
 } else if(is(x, "data.table")) {
    res <- .coercetrackingFLQuant(copy(x))
  } else {
    stop("x must be an FLmse, FLmses, or data.table")
  }
  return(res)
}

.coercetrackingFLQuant <- function(x) {
  
  # CHECK x is data.table with correct columns
  if(!is(x, "data.table") | !all(c("metric", "year", "iter", "data") %in%
    colnames(x)))
    stop("x must be a data.table with columns: (biol), metric, year, iter, data")

  # SET biol to unique if not present
  x[, biol:=ifelse(biol == "", "unique", biol)] 

  # CHANGE biol colname to unit
  setnames(x, "biol", "unit")

  # COERCE to FLQuant
  res <- as.FLQuant(x[, .(metric, year, iter, data)])

  return(res)
}

# decisions {{{

#' @examples
#' data(plesim)

decisions <- function(x, years=dimnames(tracking(x))$year, iter=NULL) {

  # EXTRACT tracking and args
  trac <- tracking(x)
  args <- args(x)

  # USE year as numeric
  years <- as.numeric(years)

  # SET iters if not given
  if(is.null(iter))
    iter <- seq(dims(x)$iter)

  # FUNCTION to compute table along years
  .table <- function(d) {

    its <- dims(d)$iter
    dmns <- dimnames(d)

    if(its == 1) {
      data.frame(metric=dmns$metric, year=dmns$year, value=prettyNum(d))
    } else {
      data.frame(metric=dmns$metric, year=dmns$year,
        value=sprintf("%s (%s)", 
          prettyNum(apply(d, 1:5, median, na.rm=TRUE)),
          prettyNum(apply(d, 1:5, mad, na.rm=TRUE))))
    }
  }

  # COMPUTE tables
  res <- lapply(years, function(y) {

    # GET advice, data and management years
    ay  <-  y
    dy <- ay - args$data_lag
    my  <- ay + args$management_lag

    # SET metrics to extract

    # data
    dmet <- c("SB.om", "SB.obs", "SB.est", "met.hcr")
    dmet <- dmet[dmet %in% dimnames(trac)$metric]

    # advice
    amet <- c("decision.hcr", "fbar.hcr", "hcr", "fbar.isys", "isys",
      "fwd", "C.om")

    amet <- amet[amet %in% dimnames(trac)$metric]

    # SUBSET metrics from tracking
    dout <- trac[dmet, ac(dy),,,, iter]
    aout <- trac[amet, ac(ay),,,, iter]

    # COMPUTE diff metrics
    mout <- trac["SB.om", ac(my),,,,iter] / trac["SB.om", ac(ay),,,,iter]
    dimnames(mout)$metric <- "diff(SB.om)"

    # BIND into single table
    rbind(.table(dout), .table(aout), .table(mout))      
  })

  do.call(cbind, res)
}
# }}}

