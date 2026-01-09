# tracking.R - DESC
# /home/mosqu003/Active/mse_FLR/msemodules/R/tracking.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# medmad {{{

medmad <- function(x) paste0(format(median(x), digits=3), " (",
  format(mad(x), digits=3) , ")")

# }}}

# inspect {{{

#' Inspect MSE Tracking Data
#'
#' The `inspect` function is designed to extract and format the tracking data from an `FLmse` object, performing aggregations and subsetting based on specified metrics and a summary function.
#'
#' @param tab An `FLmse` object or a `data.table` containing tracking data. If an `FLmse` object is supplied, the `tracking` slot is extracted.
#' @param metrics A character vector specifying the metrics to include in the output. If a single character value is given and it does not match excatly any of the contained metrics, it is used to subset using the datatable::%ilike% function. The special value `"decisions"` can be used to select all metrics from `"hcr"` onward. If `NULL`, the default, all metrics are returned.
#' @param summary A function (such as `mean` or `median`) to summarize the `data` column within the tracking `data.table` across rthe `iter` dimension. Defaults to `medmad` which returns a string with  `"Median (Median Absolute Deviation)"`.
#'
#' @details
#' The function processes tracking data in the following ways:
#'   - Extracts the tracking data from an `FLmse` object.
#'   - Arranges metrics based on the order in which they are produced inside `mp()`.
#'   - Subsets metrics based on user input or predefined criteria like `"decisions"`.
#'   - Aggregates the `data` column using the `summary` function provided.
#'   - Reshapes the processed `data.table` into a year-by-metrics format for easy reading.
#' }
#'
#' The predefined metric order is: `"om"`, `"obs"`, `"est"`, `"ind"`, `"phcr"`, `"hcr"`, `"isys"`, `"tm"`, `"iem"`, `"fb"`, `"fwd"`. Extra tracks produced inside any module are placed following the one named after the module that produced it.
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
  metord <- unique(c("year", unlist(lapply(ord, function(x)
    allmetrics[allmetrics %like% x]))))

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

  # SET columns to return
  cols <- c("year", metord)

  # RESHAPE and REORDER
  res <- dcast(tab[metric %in% metord, .(data=summary(data)), by=.(year, metric)],
    year ~ metric, value.var='data')[, ..cols]

  return(res)
}

# }}}
