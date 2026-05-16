# performance.R - DESC
# msemodules/msemodules/R/performance.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


globalVariables(c("unit"))

# CRUD: write, read, update, delete

# writePerformance {{{

#' Write performance table to file
#'
#' @param dat data.table with performance statistics.
#' @param file file to write to, defaults to 'model/performance.dat.gz' a compressed text file.
#' @return file name, invisibly.
#' @author Iago Mosqueira, WMR
#' @keywords utilities

writePerformance <- function(dat, file="model/performance.dat.gz", overwrite=FALSE) {

  # HACK to avoid method, for now
  if(is(dat, 'FLmse') | is(dat, 'FLmses')) {
    dat <- performance(dat)
  }

  # SET correct column types
  dat[, (colnames(dat)) := lapply(.SD, as.character), .SDcols = colnames(dat)]
  dat[, (c("year", "data")) := lapply(.SD, as.numeric), .SDcols = c("year", "data")]

  # ADD empty type and run if missing
  if(all(!c("type", "run") %in% names(dat))) {
    dat[, `:=`(type=character(1), run=character(1), mp=character(1)), ] 

  # SET mp from om, type and run
  } else if (is.null(dat[["mp"]]) & all(c("type", "run") %in% names(dat))) {
    dat[, mp := paste(om, type, run, sep="_")]
  }

  # SET label
  if(!"label" %in% colnames(dat)) {
    dat[, label := ifelse(mp == character(1), om, mp)]
  }

  # SET column order
  setcolorder(dat, neworder=c('om', 'type', 'run', 'mp', 'biol', 'statistic',
    'name', 'desc', 'year', 'iter', 'data'))

  # CREATE
  if(!file.exists(file) | overwrite) {

    fwrite(dat, file=file)

    invisible(TRUE)

  # ADD by substituting
  } else {

    # CHECK dat exists in file
    db <- readPerformance(file)

    # RUN anti-join on biol, statistic, year, iter, om, type & run
    db <- db[!dat, on=.(biol, statistic, year, iter, om, type, run)]

    # ADD new rows
    db <- rbind(db, dat)

    # WRITE to file
    fwrite(db, file=file)

    invisible(TRUE)
  }
}
# }}}

# readPerformance {{{

#' Read a performance table from file
#'
#' Reads a performance statistics table previously written by
#' [writePerformance()], restores column classes, keying, and column order, and
#' returns the result as a `data.table`.
#'
#' This function is the standard reader for performance tables stored on disk,
#' typically in compressed text format such as `model/performance.dat.gz`.
#'
#' @param file A character string giving the path to the performance table file.
#'   Defaults to `"model/performance.dat.gz"`.
#'
#' @return
#' A `data.table` containing performance statistics. The returned table:
#' \itemize{
#'   \item has columns coerced to the expected classes;
#'   \item is keyed by `om`, `type`, `run`, `biol`, `mp`, `statistic`, and `year`;
#'   \item has columns ordered as
#'   `om`, `type`, `run`, `mp`, `biol`, `statistic`, `name`, `desc`, `year`,
#'   `iter`, `data`;
#'   \item converts grouping columns such as `om`, `type`, `run`, `mp`, `biol`,
#'   `statistic`, and `label` to factors when present in the file.
#' }
#'
#' @details
#' The function uses [data.table::fread()] to import the table and explicitly
#' sets classes for commonly used columns in performance output:
#' \itemize{
#'   \item `type`, `run`, `mp`, and `biol` are read as character;
#'   \item `year` and `data` are read as numeric;
#'   \item `iter` is read as character.
#' }
#'
#' After reading, a key is assigned with [data.table::setkey()] to facilitate
#' fast joins and subsetting. The function also standardizes column order so
#' downstream utilities can assume a stable structure.
#'
#' If a `label` column is present in the file, it is converted to a factor along
#' with the other categorical columns.
#'
#' @seealso [mse::performance()] [writePerformance()], [summaryPerformance()], [labelPerformance()]
#'
#' @examples
#' \dontrun{
#' # Read the default performance table
#' dat <- readPerformance()
#'
#' # Read from a specific file
#' dat <- readPerformance("results/performance.dat.gz")
#'
#' # Inspect the structure
#' str(dat)
#' key(dat)
#' }
#'
#' @author Iago Mosqueira, WMR
#' @keywords file

readPerformance <- function(file="model/performance.dat.gz") {

  # READ file
  dat <- fread(file, colClasses=c(type='character', run='character',
    mp='character', biol='character', year='numeric', iter='character',
    data='numeric'))

  # SET key
  setkey(dat, om, type, run, biol, mp, statistic, year)

  # SET column order
  setcolorder(dat, neworder=c('om', 'type', 'run', 'mp', 'biol', 'statistic',
    'name', 'desc', 'year', 'iter', 'data'))

  # SET as factor
  cols <- c('om', 'type', 'run', 'mp', 'biol', 'statistic', 'label')
  dat[, (cols) := lapply(.SD, factor), .SDcols = cols]

  # RETURN
  return(dat[])
}

# }}}

# summaryPerformance {{{

#' Summarize a performance table
#'
#' Produces a compact summary of a performance statistics table, either from a
#' file on disk or from an already loaded `data.table`.
#'
#' The function reports the number of operating models (`om`), management
#' procedure types (`type`), and management procedures (`mp`), and computes a
#' grouped summary by `om`, `type`, and `run`.
#'
#' @param file Either:
#' \itemize{
#'   \item a character string giving the path to a performance table file, or
#'   \item a `data.table` already containing performance statistics.
#' }
#' Defaults to `"model/performance.dat.gz"`.
#'
#' @return
#' Invisibly returns `TRUE`. The main effect of the function is to print a brief
#' summary to the console.
#'
#' @details
#' If `file` is not a `data.table`, it is first read with [readPerformance()].
#'
#' The grouped summary computes, for each combination of `om`, `type`, and
#' `run`:
#' \itemize{
#'   \item `years`: the minimum and maximum years pasted as a range;
#'   \item `frq`: the spacing between the first two unique years, intended as an
#'     indication of temporal frequency;
#'   \item `iter`: the number of unique iterations.
#' }
#'
#' A second summary counts the number of unique values in `om`, `type`, and
#' `mp`, and prints them as a one-line overview.
#'
#' The per-group summary object is created internally but not currently printed,
#' except for the one-line counts written with [base::cat()]. The commented code
#' in the function suggests this output may later be extended to print a tree or
#' tabular summary.
#'
#' @section Current behavior and caveats:
#' \itemize{
#'   \item The function is primarily intended for quick console inspection.
#'   \item The `frq` calculation assumes at least two unique years are present
#'     within each group; otherwise it may return `numeric(0)` or `NA`.
#'   \item The grouped summary is not currently returned, only computed
#'     internally.
#' }
#'
#' @seealso [readPerformance()], [writePerformance()]
#'
#' @examples
#' \dontrun{
#' # Summarize from default file
#' summaryPerformance()
#'
#' # Summarize from a specific file
#' summaryPerformance("results/performance.dat.gz")
#'
#' # Summarize an object already in memory
#' dat <- readPerformance("results/performance.dat.gz")
#' summaryPerformance(dat)
#' }
#'
#' @author Iago Mosqueira, WMR
#' @keywords utilities

summaryPerformance <- function(file="model/performance.dat.gz") {

  #
  if(!is(file, "data.table"))
    file <- readPerformance(file)

  # TABLE by mp: years, (frq), iter, BY om, type
  res <- file[, .(
    # year range
    years=paste(min(year), max(year), sep="-"),
    # frequency
    frq=c(dist(sort(unique(as.numeric(year)))[1:2])),
    # no. iters
    iter=length(unique(iter))
    # DO by om, type & run
    ), by=.(om, type, run)]

  setorder(res, om, type, run)

  # GET summary row values
  summ <- file[, lapply(.SD, function(x) length(unique(x))),
    .SDcols = c("om", "type", "mp")] 

  # PRINT it
  cat(do.call(sprintf, c(list(fmt="- oms: %i, types: %i, mps: %i\n"), unlist(summ))))

  # PRINT tree or summary table
  # print(as.data.frame(res))

  invisible(TRUE)
}

# }}}

# labelPerformance {{{

#' @title labelPerformance
#' @description Creates a label column in a performance table
#' @param dat A performance statistics table, as returned by performance()
#' @param labels Labels to be inserted, as vector or data.frame/data.table
#' @return A labelled performance statistics data.table
#' @details
#' - 'numeric'
#' - vector
#' - data.frame or data.table
#' - NULL
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname performance

labelPerformance <- function(dat, labels=NULL) {

  # NO label, use mp | om
  if(is.null(labels)) {
    dat[, label:=ifelse(mp == character(1), om, mp)]
    return(dat[])

  # 'numeric', set as sequence in unique order for mp
  } else if(identical(labels, "numeric")) {
    labels <- data.table(mp=unique(dat[mp != character(1), mp]), 
      label=paste0("MP", seq(unique(dat[mp != character(1), mp]))))
  
  # LIST, convert to data.table
  } else if(is.list(labels)) {
    labels <- data.table(element=names(labels), label=unlist(labels))

  # SET as data.table JIC
  } else {
    labels <- data.table(labels)
  }

  # GET dims
  dimdat <- dim(dat)

  # CREATE tmp column to match mp | om
  dat[, element:=ifelse(mp == "", as.character(om), as.character(mp))]

  # MERGE new labels on matching element
  dat <- merge(dat[, !"label"], labels[element %in% unique(dat$element)],
    by="element", all=TRUE)

  # TODO dat <- dat[, !"label"][labels, on = .(element = element), roll = TRUE, nomatch=0]

  # SET NA to empty string
  dat[, label:=ifelse(is.na(label), element, label)]

  # DROP tmp column
  dat[, element:=NULL]
  
  # SET as factor, OM labels (no mp) first
  levs <- c(dat[mp == character(1),
    unique(label)], sort(dat[mp != character(1), unique(label)]))
 
  dat[, label := factor(label, levels=levs)]

  # CHECK dims
  if(!identical(dim(dat), dimdat))
    warning("Missmatch in dimensions of tables, check output.")

  # END
  return(dat[])
}
# }}}

# setLabelPerformance {{{

#' @title setLabelPerofrmance
#' @description FUNCTION_DESCRIPTION
#' @param file File where performance table is stored, default: 'model/performance.dat.gz'
#' @return Invisible updates the table in file
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname labelPerformance

setLabelPerformance <- function(file="model/performance.dat.gz", labels) {

  dat <- readPerformance(file)

  dat <- labelPerformance(dat, labels)

  writePerformance(dat, file=file, overwrite=TRUE)
}
# }}}

# periodsPerformance {{{

#' Aggregate performance statistics by period
#'
#' Computes mean performance values over user-defined periods and returns a
#' summarized performance table.
#'
#' This function is useful for collapsing annual performance time series into
#' broader periods, such as short-, medium-, and long-term intervals, while
#' preserving the main grouping structure of the original performance table.
#'
#' @param x A `data.table` containing performance statistics. It must include at
#'   least the columns `year`, `data`, `type`, `mp`, `statistic`, `name`,
#'   `desc`, and `iter`. If a `label` column is present, it is retained in the
#'   grouping structure.
#' @param periods A vector or list defining the periods over which means should
#'   be computed. Each element should contain the years belonging to one period.
#'   Named elements are used as period labels; unnamed elements are labelled
#'   automatically from the year range.
#'
#' @return
#' A `data.table` with one row per grouping combination and period, containing:
#' \itemize{
#'   \item the grouping variables inherited from `x`;
#'   \item `period`, the name of the period;
#'   \item `year`, a compact character representation of the year span;
#'   \item `data`, the mean of `data` over the years in that period.
#' }
#'
#' @details
#' The `periods` argument is coerced to a list. For each period:
#' \itemize{
#'   \item if the element has length greater than one, the display label in the
#'     output `year` column is built as `"start-end"`, where `end` uses the last
#'     two digits of the final year;
#'   \item if the element has length one, that year is used directly.
#' }
#'
#' If period names are missing, they are replaced by the derived year labels.
#'
#' Means are computed with `na.rm=TRUE`, so missing values in `data` are ignored.
#'
#' When `x` contains a `label` column, aggregation is done by:
#' `type`, `mp`, `label`, `statistic`, `name`, `desc`, and `iter`.
#'
#' Otherwise aggregation is done by:
#' `type`, `mp`, `statistic`, `name`, `desc`, and `iter`.
#'
#' @section Output structure:
#' The resulting table is intended to remain compatible with downstream plotting
#' or reporting workflows based on summarized performance indicators.
#'
#' @seealso [labelPerformance()], [summaryPerformance()]
#'
#' @examples
#' \dontrun{
#' dat <- readPerformance()
#'
#' # Define named periods
#' per <- list(
#'   short = 2021:2025,
#'   medium = 2026:2035,
#'   long = 2036:2050
#' )
#'
#' res <- periodsPerformance(dat, per)
#'
#' # Unnamed periods will be labelled automatically
#' res2 <- periodsPerformance(dat, list(2021:2025, 2026:2030))
#' }
#'
#' @author Iago Mosqueira, WMR
#' @keywords manip
#' @export

periodsPerformance <- function(x, periods) {

  # COERCE to list
  periods <- as.list(periods)
 
  years <- unlist(lapply(periods, function(x) {
    if(length(x) > 1)
      paste(x[1], substr(rev(x)[1], 3, 4), sep="-")
    else
      x
  }))

  # ASSIGN names if missing
  names(periods)[names(periods) == character(1)] <-
    years[names(periods) == character(1)]

  # COMPUTE means per period by label or mp
  if("label" %in% colnames(x)) {
    res <- rbindlist(Map(function(pe, na, ye) {
      x[year %in% pe, .(data=mean(data, na.rm=TRUE), period=na, year=ye),
      by=.(type, mp, label, statistic, name, desc, iter)]},
      pe=periods, na=names(periods), ye=years))
  } else {
    res <- rbindlist(Map(function(pe, na, ye) {
      x[year %in% pe, .(data=mean(data, na.rm=TRUE), period=na, year=ye),
      by=.(type, mp, statistic, name, desc, iter)]},
      pe=periods, na=names(periods), ye=years))
  }

  return(res)
}
# }}}

# extractPerformance {{{

#' Extract performance data for one or more management procedures
#'
#' Subsets a performance table to return the time series corresponding to one or
#' more management procedures (`mp`) together with the associated historical
#' operating model (`om`) rows.
#'
#' This is useful when comparing the performance of a selected management
#' procedure against its corresponding operating model baseline.
#'
#' @param dat A performance statistics `data.table` containing at least the
#'   columns `om` and `mp`.
#' @param mp A character string used to match management procedures. Matching is
#'   performed with `%like%`, so partial strings or patterns can be used.
#'
#' @return
#' A subset of `dat` containing:
#' \itemize{
#'   \item rows for management procedures whose `mp` matches the supplied
#'     pattern; and
#'   \item rows for the corresponding operating model(s), identified as rows
#'     where `om` matches the selected subset and `mp` is empty (`""`).
#' }
#'
#' @details
#' The function first identifies all rows whose `mp` matches the supplied
#' pattern using [data.table::like()] syntax via `%like%`. From this subset, it
#' extracts the unique matching management procedures and associated operating
#' models. It then returns all rows in the original table that belong either to:
#' \itemize{
#'   \item one of the matched management procedures, or
#'   \item the corresponding operating model baseline rows.
#' }
#'
#' This design makes it easy to compare projected management procedure
#' trajectories against the historical or reference operating model records.
#'
#' @section Current limitations:
#' \itemize{
#'   \item The function contains a TODO note indicating future support for
#'     parsing multiple management procedure inputs more explicitly.
#'   \item Matching is pattern-based rather than exact, so care should be taken
#'     when `mp` names overlap.
#'   \item The function assumes baseline operating model rows are identified by
#'     an empty string in the `mp` column.
#' }
#'
#' @seealso [labelPerformance()], [periodsPerformance()]
#'
#' @examples
#' \dontrun{
#' dat <- readPerformance()
#'
#' # Extract rows for one MP pattern and its corresponding OM
#' sub <- extractPerformance(dat, "HCR1")
#'
#' # Because matching uses %like%, partial patterns can be used
#' sub2 <- extractPerformance(dat, "MP")
#' }
#'
#' @author Iago Mosqueira, WMR
#' @keywords manip
#' @export

extractPerformance <- function(dat, mp) {

  # TODO: PARSE multiple MPs and match each OM

  # ASSIGN to avoid column match 
  smp <- mp

  dat[mp %like% smp]

  # FIND mps & om
  sub <- dat[mp %like% smp]
  mps <- sub[, as.character(unique(mp))]
  oms <- sub[, as.character(unique(om))]

  # RETURN subset, om + mp
  return(dat[om %in% oms & mp %in% c("", mps)])
}
# }}}

# getPerformance {{{

#' Build a performance table from operating model files
#'
#' Reads a set of serialized operating model objects from disk, computes
#' performance statistics for each, and combines the results into a single
#' `data.table`.
#'
#' @param path A character string giving the directory containing serialized
#'   operating model files.
#' @param pattern A character string giving the file pattern used by
#'   [base::list.files()] to identify files to read. Defaults to `"*.rds"`.
#' @param fy The final year to retain when windowing each operating model before
#'   computing performance. Passed to [window()].
#' @param ... Additional arguments passed to [performance()].
#'
#' @return
#' A `data.table` containing the row-bound performance statistics for all files
#' found in `path` matching `pattern`. An `om` column is added using the file
#' name with the `.rds` extension removed.
#'
#' @details
#' For each file found:
#' \enumerate{
#'   \item the object is read using [readRDS()];
#'   \item the `$om` component is extracted;
#'   \item the object is trimmed to `fy` using [window()];
#'   \item [performance()] is called on the resulting object;
#'   \item an `om` column is added from the file basename.
#' }
#'
#' The results are combined using [data.table::rbindlist()].
#'
#' Warnings generated while computing performance are suppressed using
#' [base::suppressWarnings()], which can be useful in batch processing but may
#' hide potentially important diagnostics.
#'
#' @section Assumptions:
#' \itemize{
#'   \item Each `.rds` file contains an object with an `$om` component.
#'   \item That component is compatible with [window()] and [performance()].
#'   \item File names uniquely identify the operating model and can be used as
#'     the `om` label.
#' }
#'
#' @seealso [getMSEPerformance()], [writePerformance()], [performance()]
#'
#' @examples
#' \dontrun{
#' # Build a performance table from OM files in a directory
#' dat <- getOMPerformance("om", fy = 2025)
#'
#' # Pass additional arguments on to performance()
#' dat <- getOMPerformance("om", fy = 2025, probs = c(0.1, 0.5, 0.9))
#' }
#'
#' @author Iago Mosqueira, WMR
#' @keywords file
#' @export

getOMPerformance <- function(path, pattern="*.rds", fy, ...) {
  return(rbindlist(lapply(list.files(path, pattern, full.names=TRUE), function(i)
    suppressWarnings(performance(window(readRDS(i)$om, end=fy), ...)[,
      om:=sub('.rds', '', basename(i))])[])))
}

#' Build a performance table from MSE result files
#'
#' Reads a set of serialized MSE result files from disk and combines their
#' performance statistics into a single `data.table`.
#'
#' @param path A character string giving the directory containing serialized MSE
#'   result files.
#' @param pattern A character string giving the file pattern used by
#'   [base::list.files()] to identify files to read. Defaults to `"*.rds"`.
#'
#' @return
#' A `data.table` created by row-binding performance results from all matching
#' files.
#'
#' @details
#' For each file found in `path` matching `pattern`, the object is read using
#' [readRDS()]. Then:
#' \itemize{
#'   \item if the object already inherits from `data.table`, it is returned as-is;
#'   \item otherwise, [performance()] is called on the object to derive the
#'     performance table.
#' }
#'
#' All resulting tables are combined with [data.table::rbindlist()] using
#' `fill=TRUE`, so files with slightly different column structures can still be
#' merged.
#'
#' This function is intended as a convenient batch importer for previously saved
#' MSE outputs or already extracted performance tables.
#'
#' @section Typical use cases:
#' \itemize{
#'   \item combining performance output from multiple management procedures;
#'   \item importing a mixture of saved raw MSE objects and precomputed
#'     `data.table` performance summaries;
#'   \item preparing inputs for [writePerformance()] or downstream summaries.
#' }
#'
#' @seealso [getOMPerformance()], [readPerformance()], [writePerformance()]
#'
#' @examples
#' \dontrun{
#' # Read and combine performance tables from a folder of MSE results
#' dat <- getMSEPerformance("mse")
#'
#' # Use a custom file pattern
#' dat <- getMSEPerformance("mse", pattern = "scenario_.*\\.rds$")
#' }
#'
#' @author Iago Mosqueira, WMR
#' @keywords file
#' @export

getMSEPerformance <- function(path, pattern="*.rds") {
  return(rbindlist(lapply(list.files(path, pattern, full.names=TRUE), function(i) {
    dat <- readRDS(i)
    if(is(dat, "data.table"))
      return(dat)
    else 
      return(performance(dat))
    }
  ), fill=TRUE))
}
# }}}

# performanceFLQuants {{{

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
#'   quant <- performanceFLQuant(perf_dt)
#'   quants <- performanceFLQuants(perf_dt)
#' }
#'
#' @seealso [mse::performance()], [mse::FLmse-class], [mse::FLmses-class]
#'
#' @keywords manip

performanceFLQuant <- function(x) {
  
  res <- performanceFLQuants(x)

  if(is.list(res) && length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}

#' @rdname performanceFLQuant

performanceFLQuants <- function(x) {

  # USE poor mans' dispatch
  if(is(x, "FLmse") | is(x, "FLmses")) {
    x <- performance(x)
  }

  # CHECK x is data.table with correct columns
  if(!is(x, "data.table") | !all(c("statistic", "year", "iter", "data", "mp") %in%
    colnames(x)))
    stop("x must be a data.table with columns: statistic, year, iter, data, mp")
  
  # CHECK no. of biols
  if("biol" %in% colnames(x)) {
    setnames(x, "biol", "unit")
   
    # BUILD FLQuant per mp, biol as unit
    res <- lapply(split(x[, .(statistic, year, unit, iter, data, mp)], by="mp"),
      function(x) {
        as.FLQuant(x[, .(statistic, year, unit, iter, data)])
      }
    )
  } else {
    # BUILD FLQuant per mp
    res <- lapply(split(x[, .(statistic, year, iter, data, mp)], by="mp"),
      function(x) {
        as.FLQuant(x[, .(statistic, year, iter, data)])
      }
    )
  }

  return(res)
}
# }}}
