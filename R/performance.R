# performance.R - DESC
# msemodules/msemodules/R/performance.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


globalVariables(c("unit"))

# .atomicWrite {{{

.atomicWrite <- function(db, file) {
  tmp <- paste0(file, ".tmp")
  on.exit(unlink(tmp), add = TRUE)
  fwrite(db, file = tmp)
  file.rename(tmp, file)
  invisible(TRUE)
}
# }}}

# validatePerformance {{{

#' @describeIn writePerformance Validate a performance statistics table,
#'   checking for required columns, correct numeric types, NA values in key
#'   columns, and duplicate rows on the composite primary key.
#' @keywords utilities

validatePerformance <- function(dat) {

  required <- c("om", "statistic", "year", "iter", "data")
  missing  <- setdiff(required, names(dat))
  if(length(missing))
    stop("Missing required columns: ", paste(missing, collapse=", "))

  if(!is.numeric(dat[["year"]]))
    stop("Column 'year' must be numeric.")
  if(!is.numeric(dat[["data"]]))
    stop("Column 'data' must be numeric.")

  key_cols <- intersect(c("om", "type", "run", "mp", "biol", "statistic",
    "year", "iter"), names(dat))
  na_counts <- dat[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = key_cols]
  bad <- names(na_counts)[unlist(na_counts) > 0]
  if(length(bad))
    warning("NA values found in key column(s): ", paste(bad, collapse=", "))

  dups <- sum(duplicated(dat[, ..key_cols]))
  if(dups > 0)
    warning(dups, " duplicate row(s) on primary key - table may not be unique after write.")

  invisible(TRUE)
}
# }}}

# hasPerformance {{{

#' @describeIn writePerformance Test whether a stored performance table
#'   contains rows matching the supplied identifiers, without reading the full
#'   table into memory.
#' @param om,type,run Character scalar identifiers used to filter the stored
#'   table. Any combination may be supplied; unspecified arguments are ignored.
#' @keywords file

hasPerformance <- function(file="model/performance.dat.gz", om=NULL, type=NULL, run=NULL) {

  if(!file.exists(file))
    return(FALSE)

  db <- fread(file, select=c("om", "type", "run"))

  # COPY args to unambiguous locals — avoids column/arg name clash in [...]
  .om   <- om
  .type <- type
  .run  <- run

  if(!is.null(.om))   db <- db[om   == .om]
  if(!is.null(.type)) db <- db[type == .type]
  if(!is.null(.run))  db <- db[run  == .run]

  return(nrow(db) > 0L)
}
# }}}

# listPerformance {{{

#' @describeIn writePerformance Return a compact catalogue of the
#'   \code{(om, type, run, mp)} combinations stored in a performance file,
#'   together with row count, year range, and number of distinct iterations
#'   without loading the \code{data} column.
#' @keywords file

listPerformance <- function(file="model/performance.dat.gz") {

  if(!file.exists(file))
    stop("File not found: ", file)

  db <- fread(file, select=c("om", "type", "run", "mp", "year", "iter"))

  res <- db[, .(
    rows  = .N,
    years = paste(min(year), max(year), sep="-"),
    iters = length(unique(iter))
  ), by=.(om, type, run, mp)]

  setorder(res, om, type, run, mp)

  return(res[])
}
# }}}

# deletePerformance {{{

#' @describeIn writePerformance Delete rows from a stored performance table
#'   that match the supplied identifiers.  When \code{dry_run = TRUE} the rows
#'   that \emph{would} be deleted are returned without modifying the file.
#' @param dry_run Logical. If \code{TRUE} (default \code{FALSE}), return the
#'   rows that would be deleted without writing any changes.
#' @keywords file

deletePerformance <- function(file="model/performance.dat.gz",
  om=NULL, type=NULL, run=NULL, mp=NULL, dry_run=FALSE) {

  if(!file.exists(file))
    stop("File not found: ", file)

  db <- readPerformance(file)

  # COPY args to unambiguous locals
  .om   <- om
  .type <- type
  .run  <- run
  .mp   <- mp

  mask <- rep(TRUE, nrow(db))
  if(!is.null(.om))   mask <- mask & (as.character(db[["om"]])   %in% .om)
  if(!is.null(.type)) mask <- mask & (as.character(db[["type"]]) %in% .type)
  if(!is.null(.run))  mask <- mask & (as.character(db[["run"]])  %in% .run)
  if(!is.null(.mp))   mask <- mask & (as.character(db[["mp"]])   %in% .mp)

  deleted <- db[mask]
  message(sum(mask), " row(s) selected for deletion.")

  if(dry_run)
    return(deleted[])

  .atomicWrite(db[!mask], file)
  invisible(deleted[])
}
# }}}

# diffPerformance {{{

#' @describeIn writePerformance Compare an in-memory performance table against
#'   a stored file and categorise every row as \emph{new}, \emph{replace}
#'   (key exists in file), or \emph{unchanged} (key exists and values are
#'   identical).  Returns a named list of three \code{data.table}s; the file
#'   is not modified.
#' @keywords file

diffPerformance <- function(dat, file="model/performance.dat.gz") {

  if(!file.exists(file)) {
    message("File does not exist, all rows are new.")
    return(list(new=dat[], replace=dat[0], unchanged=dat[0]))
  }

  key_cols <- intersect(c("om", "type", "run", "biol", "statistic", "year", "iter"),
    names(dat))

  db <- readPerformance(file)

  # NEW rows: key not in db
  new_rows <- dat[!db, on=key_cols]

  # Rows whose key IS in db
  existing <- dat[db, on=key_cols, nomatch=NULL]

  # UNCHANGED: key + data match
  unchanged <- existing[db, on=c(key_cols, "data"), nomatch=NULL]

  # REPLACE: key matches but data differs
  replace <- existing[!unchanged, on=key_cols]

  message(nrow(new_rows), " new, ", nrow(replace), " replace, ",
    nrow(unchanged), " unchanged.")

  return(list(new=new_rows[], replace=replace[], unchanged=unchanged[]))
}
# }}}

# writePerformance {{{

#' Read and write performance statistics tables
#'
#' A family of utilities for persisting, retrieving, labelling, summarizing,
#' and reshaping performance statistics tables produced by
#' [mse::performance()].
#'
#' \describe{
#'   \item{[validatePerformance()]}{Checks required columns, numeric types,
#'     NA values in key columns, and key uniqueness; stops or warns on failure.}
#'   \item{[hasPerformance()]}{Returns \code{TRUE} if the file contains rows
#'     matching the supplied identifiers, without loading the full table.}
#'   \item{[listPerformance()]}{Returns a catalogue \code{data.table} of
#'     \code{(om, type, run, mp)} combinations with row counts and year ranges.}
#'   \item{[deletePerformance()]}{Deletes rows matching supplied identifiers;
#'     supports a \code{dry_run} preview mode.}
#'   \item{[diffPerformance()]}{Compares an in-memory table against the stored
#'     file and categorises rows as new, replace, or unchanged.}
#'   \item{[writePerformance()]}{Serializes a table to disk, merging with any
#'     existing file.}
#'   \item{[readPerformance()]}{Reads a stored table and restores column
#'     classes, key, and factor levels.}
#'   \item{[summaryPerformance()]}{Prints a compact console summary.}
#'   \item{[labelPerformance()]}{Adds or replaces the \code{label} column.}
#'   \item{[setLabelPerformance()]}{Reads, re-labels, and overwrites a stored
#'     table.}
#'   \item{[periodsPerformance()]}{Aggregates annual statistics into
#'     user-defined periods.}
#'   \item{[extractPerformance()]}{Subsets to selected MPs together with their
#'     OM baseline rows.}
#'   \item{[getOMPerformance()]}{Batch-reads OM files and returns a combined
#'     performance table.}
#'   \item{[getMSEPerformance()]}{Batch-reads MSE result files and returns a
#'     combined performance table.}
#' }
#'
#' @param dat A `data.table` containing performance statistics, as returned by
#'   [mse::performance()], or an object of class `FLmse` or `FLmses` from
#'   which performance statistics can be extracted automatically.
#' @param file A character string giving the path to the performance table
#'   file. Defaults to `"model/performance.dat.gz"`. The `.gz` extension
#'   causes [data.table::fwrite()] / [data.table::fread()] to use gzip
#'   compression transparently.
#' @param overwrite Logical. If `TRUE`, any existing file at `file` is
#'   overwritten entirely. If `FALSE` (default), existing rows not conflicting
#'   with `dat` are preserved and new rows are appended; conflicts are resolved
#'   by an anti-join on `biol`, `statistic`, `year`, `iter`, `om`, `type`,
#'   and `run`.
#' @param labels The labelling specification used by [labelPerformance()] and
#'   [setLabelPerformance()]. Accepted forms:
#'   \describe{
#'     \item{`NULL` (default)}{Each row is labelled by its `mp` value if one
#'       is present, or by `om` otherwise.}
#'     \item{`"numeric"`}{Management procedures are assigned sequential labels
#'       `"MP1"`, `"MP2"`, ... in the order they appear.}
#'     \item{A named `list`}{Names are matched against `mp` (or `om`) and the
#'       values are used as labels.}
#'     \item{A `data.frame` or `data.table`}{Must contain columns `element`
#'       and `label`; `element` is matched against `mp` (or `om`).}
#'   }
#' @param x A `data.table` containing performance statistics for
#'   [periodsPerformance()], as returned by [mse::performance()] or
#'   [readPerformance()]. Must contain at least `year`, `data`, `type`, `mp`,
#'   `statistic`, `name`, `desc`, and `iter`.
#' @param periods A vector or list whose elements define the years belonging to
#'   each period (used by [periodsPerformance()]). Named elements use the name
#'   as the period label; unnamed elements are labelled automatically from the
#'   year range (e.g. `"2026-35"`).
#' @param mp A character value (or vector) used by [extractPerformance()].
#'   Length-1 values are matched with [data.table::like()] (`%like%`) so
#'   partial strings and regular-expression patterns are accepted; vectors are
#'   matched exactly with `%in%`.
#' @param path A character string giving the directory containing serialized
#'   files, used by [getOMPerformance()] and [getMSEPerformance()].
#' @param pattern A character string passed to [base::list.files()] to filter
#'   which files are read. Defaults to `"*.rds"`.
#' @param fy The final year to which each operating model is trimmed with
#'   [FLCore::window()] before performance is computed
#'   ([getOMPerformance()] only).
#' @param ... Additional arguments forwarded to [mse::performance()]
#'   ([getOMPerformance()] only).
#'
#' @return
#' \describe{
#'   \item{[writePerformance()], [setLabelPerformance()]}{Invisibly `TRUE`.}
#'   \item{[deletePerformance()]}{Invisibly, a \code{data.table} of the deleted
#'     rows. When \code{dry_run = TRUE}, the rows that would be deleted are
#'     returned visibly without modifying the file.}
#'   \item{[diffPerformance()]}{A named list with elements \code{new},
#'     \code{replace}, and \code{unchanged}, each a \code{data.table}.}
#'   \item{[readPerformance()]}{A `data.table` keyed by `om`, `type`, `run`,
#'     `biol`, `mp`, `statistic`, and `year`, with grouping columns as
#'     factors and column order fixed as `om`, `type`, `run`, `mp`, `biol`,
#'     `statistic`, `name`, `desc`, `year`, `iter`, `data`.}
#'   \item{[summaryPerformance()]}{Invisibly `TRUE` (called for its printed
#'     side-effect).}
#'   \item{[labelPerformance()]}{A copy of `dat` with a `label` factor column
#'     added or replaced; OM-only rows appear first in the factor levels.}
#'   \item{[periodsPerformance()]}{A `data.table` with one row per grouping
#'     combination and period, containing `period`, `year` (compact range
#'     string), `data` (period mean), and the inherited grouping columns.}
#'   \item{[extractPerformance()]}{A subset of `dat` containing the matched MP
#'     rows and their corresponding OM baseline rows (`mp == ""`).}
#'   \item{[getOMPerformance()], [getMSEPerformance()]}{A `data.table` produced
#'     by row-binding the performance results from all matching files.}
#' }
#'
#' @details
#' **`writePerformance`**
#'
#' Before writing, all columns are coerced to character and then `year` and
#' `data` are restored to numeric to ensure a consistent serialization format.
#' Optional columns are normalized: if neither `type` nor `run` is present,
#' empty character columns for `type`, `run`, and `mp` are added; if `mp` is
#' absent but `type` and `run` are present, `mp` is constructed by pasting
#' `om`, `type`, and `run` with `"_"`; if `label` is absent it is set to `mp`
#' when an MP is present, otherwise to `om`. Column order is standardized to
#' `om`, `type`, `run`, `mp`, `biol`, `statistic`, `name`, `desc`, `year`,
#' `iter`, `data` before writing.
#'
#' **`readPerformance`**
#'
#' Reading is performed with [data.table::fread()] with explicit `colClasses`
#' to avoid ambiguous type inference (e.g. when `iter` contains non-numeric
#' labels or `mp` is an empty string). The key set by [data.table::setkey()]
#' enables efficient subsetting and is assumed by several downstream functions.
#' The `label` column is optional; when absent it is not created, use
#' [labelPerformance()] to add it after reading.
#'
#' **`periodsPerformance`**
#'
#' `periods` is coerced to a list internally. The compact year label uses only
#' the last two digits of the final year (`2026:2035` to `"2026-35"`). Missing
#' values in `data` are silently ignored (`na.rm = TRUE`). When `x` contains a
#' `label` column the grouping includes it; otherwise grouping is by `type`,
#' `mp`, `statistic`, `name`, `desc`, and `iter`.
#'
#' **`extractPerformance`**
#'
#' Operating model baseline rows are identified by an empty string in `mp`.
#' For scalar `mp`, matching uses `%like%`; use anchored patterns (e.g.
#' `"^HCR1$"`) when MP names share substrings. For vectors, exact matching
#' (`%in%`) is used.
#'
#' **`getOMPerformance`**
#'
#' Each `.rds` file is expected to contain a list-like object with an `$om`
#' component compatible with [FLCore::window()] and [mse::performance()].
#' Warnings from [mse::performance()] are suppressed with
#' [base::suppressWarnings()]; re-run individual files interactively if
#' unexpected results appear.
#'
#' **`getMSEPerformance`**
#'
#' If a file's deserialized object already inherits from `data.table` it is
#' returned directly, allowing pre-computed performance summaries to be mixed
#' with raw MSE objects. `fill = TRUE` in [data.table::rbindlist()] tolerates
#' minor column-structure differences across files.
#'
#' @seealso [mse::performance()], [mse::FLmse-class], [mse::FLmses-class],
#'   [FLCore::window()], [performanceFLQuants()]
#'
#' @examples
#' \dontrun{
#' ## writePerformance / readPerformance
#' writePerformance(perf_dat)
#' dat <- readPerformance()
#' dat <- readPerformance("results/performance.dat.gz")
#'
#' ## validatePerformance
#' validatePerformance(perf_dat)
#'
#' ## hasPerformance
#' hasPerformance()
#' hasPerformance(om = "ple.27.420")
#' if(!hasPerformance(om = "ple.27.420", run = "r01"))
#'   writePerformance(perf_dat)
#'
#' ## listPerformance
#' listPerformance()
#' listPerformance("results/performance.dat.gz")
#'
#' ## diffPerformance
#' d <- diffPerformance(perf_dat)
#' d$replace
#' if(nrow(d$replace) == 0) writePerformance(perf_dat)
#'
#' ## deletePerformance
#' deletePerformance(om = "ple.27.420", run = "r01", dry_run = TRUE)
#' deletePerformance(om = "ple.27.420", run = "r01")
#' deletePerformance(mp = "hcr_Fmsy")
#'
#' ## summaryPerformance
#' summaryPerformance()
#' summaryPerformance(dat)
#'
#' ## labelPerformance
#' dat <- labelPerformance(dat)
#' dat <- labelPerformance(dat, labels = "numeric")
#' dat <- labelPerformance(dat, labels = list(
#'   hcr01 = "Conservative HCR", hcr02 = "Moderate HCR"))
#'
#' ## setLabelPerformance
#' setLabelPerformance(labels = "numeric")
#' setLabelPerformance("results/performance.dat.gz",
#'   labels = list(hcr01 = "Low F", hcr02 = "High F"))
#'
#' ## periodsPerformance
#' res <- periodsPerformance(dat, periods = list(
#'   short  = 2021:2025,
#'   medium = 2026:2035,
#'   long   = 2036:2050))
#'
#' ## extractPerformance
#' sub <- extractPerformance(dat, "HCR1")
#' sub <- extractPerformance(dat, "^HCR1$")
#'
#' ## getOMPerformance
#' dat <- getOMPerformance("om", fy = 2025)
#' dat <- getOMPerformance("om", fy = 2025, probs = c(0.1, 0.5, 0.9))
#'
#' ## getMSEPerformance
#' dat <- getMSEPerformance("mse")
#' dat <- getMSEPerformance("mse", pattern = "scenario_.*\\.rds$")
#' }
#'
#' @author Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#' @name writePerformance
#' @keywords file

writePerformance <- function(dat, file="model/performance.dat.gz", overwrite=FALSE) {

  # COPY to avoid direct changes
  dat <- copy(dat)

  # HACK to avoid method, for now
  if(is(dat, 'FLmse') | is(dat, 'FLmses')) {
    dat <- performance(dat)
  }

  # ADD empty biol if missing (single-stock / FLom case)
  if(!"biol" %in% names(dat))
    dat[, biol := ""]

  # VALIDATE table
  validatePerformance(dat)

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

    .atomicWrite(dat, file)

    invisible(TRUE)

  # ADD by substituting
  } else {

    # CHECK dat exists in file
    db <- readPerformance(file)

    # RUN anti-join on biol, statistic, year, iter, om, type & run
    db <- db[!dat, on=.(biol, statistic, year, iter, om, type, run)]

    # ADD new rows
    new_cols <- setdiff(names(dat), names(db))
    old_cols <- setdiff(names(db),  names(dat))

    # WARN of new and old columns
    if(length(new_cols))
      message("New columns in dat not in file (filled with NA): ", paste(new_cols, collapse=", "))

    if(length(old_cols))
      message("Columns in file not in dat (filled with NA): ",    paste(old_cols, collapse=", "))
    
    # BIND tables
    db <- data.table::rbindlist(list(db, dat), fill = TRUE)

    # WRITE to file
    .atomicWrite(db, file)

    invisible(TRUE)
  }
}
# }}}

# readPerformance {{{

#' @describeIn writePerformance Read a performance statistics table from disk,
#'   restoring column classes, data.table key, column order, and factor levels.
#' @keywords file

readPerformance <- function(file="model/performance.dat.gz") {

  # FREAD forcing column classes
  dat <- fread(file, colClasses=c(type='character', run='character',
    mp='character', biol='character', year='numeric', iter='character',
    data='numeric'))

  # SET column order FIRST
  setcolorder(dat, neworder=intersect(c('om', 'type', 'run', 'mp', 'biol',
    'statistic', 'name', 'desc', 'year', 'iter', 'data'), names(dat)))

  # CONVERT grouping columns to factor
  cols <- intersect(c('om', 'type', 'run', 'mp', 'biol', 'statistic', 'label'),
    names(dat))
  dat[, (cols) := lapply(.SD, factor), .SDcols = cols]

  # SET key AFTER factor conversion — key is on numeric/character columns only
  setkey(dat, statistic, year)

  return(dat[])
}

# }}}

# summaryPerformance {{{

#' @describeIn writePerformance Print a compact console summary of the number
#'   of operating models, MP types, and management procedures, together with
#'   per-group year ranges and iteration counts.
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
    iter=length(unique(iter)),
    # no. statistics
    statistics=length(unique(statistic))
    # DO by om, type & run
    ), by=.(om, type, label)]

  setorder(res, om, type, run)

  # PRINT tree or summary table
  print(as.data.frame(res))

  invisible(res)
}

# }}}

# labelPerformance {{{

#' @describeIn writePerformance Add or replace the \code{label} factor column
#'   in a performance statistics table; OM-only rows appear first in the factor
#'   levels.
#' @keywords manip

labelPerformance <- function(dat, labels=NULL) {

  # COPY to avoid modifying caller's object
  dat <- copy(dat)

  # COERCE factor columns to character FIRST — must happen in all branches
  if(is.factor(dat[["mp"]])) dat[, mp := as.character(mp)]
  if(is.factor(dat[["om"]])) dat[, om := as.character(om)]
  if(is.factor(dat[["label"]])) dat[, label := as.character(label)]

  # HELPER: set label as factor with OM rows first
  .setLabelFactor <- function(dat) {
    levs <- c(
      unique(dat[mp == "", label]),
      sort(unique(dat[mp != "", label])))
    dat[, label := factor(label, levels=levs)]
    dat
  }

  # NO label, use mp | om
  if(is.null(labels)) {
    dat[, label := ifelse(mp == "", om, mp)]
    return(.setLabelFactor(dat)[])

  # 'numeric': sequential MP1, MP2, ... in order of appearance
  } else if(identical(labels, "numeric")) {
    mps <- unique(dat[mp != "", mp])
    labels <- data.table(element=mps, label=paste0("MP", seq_along(mps)))

  # LIST: convert to data.table with 'element' column
  } else if(is.list(labels)) {
    labels <- data.table(element=names(labels), label=unlist(labels))

  # data.frame/data.table: ensure correct class
  } else {
    labels <- as.data.table(labels)
  }

  # CHECK labels has required columns
  if(!all(c("element", "label") %in% names(labels)))
    stop("'labels' must have columns 'element' and 'label'.")

  # GET dims for post-merge check
  dimdat <- dim(dat)

  # CREATE tmp column: mp if set, otherwise om
  dat[, element := ifelse(mp == "", om, mp)]

  # DROP existing label column if present
  if("label" %in% names(dat))
    dat[, label := NULL]

  # MERGE new labels (left join — unmatched rows keep their element as label)
  dat <- merge(dat, labels[element %in% unique(dat$element)],
    by="element", all.x=TRUE)

  # SET NA labels to element value
  dat[, label := ifelse(is.na(label), element, label)]

  # DROP tmp column
  dat[, element := NULL]

  # CHECK dims unchanged
  if(!identical(dim(dat), dimdat))
    warning("Mismatch in dimensions of table, check output.")

  return(.setLabelFactor(dat)[])
}
# }}}

# setLabelPerformance {{{

#' @describeIn writePerformance Convenience wrapper that reads a stored
#'   performance table, updates its \code{label} column via
#'   [labelPerformance()], and writes the result back, replacing the previous
#'   contents.
#' @keywords file

setLabelPerformance <- function(file="model/performance.dat.gz", labels) {

  dat <- readPerformance(file)

  dat <- labelPerformance(dat, labels)

  writePerformance(dat, file=file, overwrite=TRUE)
}
# }}}

# periodsPerformance {{{

#' @describeIn writePerformance Collapse annual performance statistics into
#'   broader user-defined periods by computing the mean of \code{data} over
#'   the years belonging to each period.
#' @keywords manip

periodsPerformance <- function(x, periods=list(), add=FALSE, ...) {

  # COLLECT any named ... args (e.g. short=2026:2030, long=2031:2040)
  dots <- list(...)

  # MERGE with periods list — dots take precedence on name clash
  periods <- c(periods, dots)

  # MUST have at least one period
  if(length(periods) == 0)
    stop("Supply at least one period, e.g. periodsPerformance(x, short=2026:2030).")

  # COERCE to list
  periods <- as.list(periods)

  # BUILD compact year labels: 2026:2035 -> "2026-35", single year -> "2030"
  years <- unlist(lapply(periods, function(p) {
    if(length(p) > 1)
      paste(p[1], substr(as.character(rev(p)[1]), 3, 4), sep="-")
    else
      as.character(p)
  }))

  end <- unlist(lapply(periods, function(p) {
    p[length(p)]
  }))

  # ASSIGN names to unnamed periods from year label
  nms <- names(periods)
  if(is.null(nms))
    nms <- rep("", length(periods))
  nms[nms == ""] <- years[nms == ""]
  names(periods) <- nms

  # COMPUTE means per period
  grp <- if("label" %in% colnames(x))
    quote(.(om, biol, type, run, mp, label, statistic, name, desc, iter))
  else
    quote(.(om, biol, type, run, mp, statistic, name, desc, iter))

  res <- rbindlist(Map(function(pe, na, ye, en) {
    x[year %in% pe,
      .(data=mean(data, na.rm=TRUE), period=na, year=as.numeric(NA), years=ye),
      by=eval(grp)]
    }, pe=periods, na=names(periods), ye=years, en=end))

  # ADD
  if(add) {
    res <- rbindlist(list(x, res), fill=TRUE)
  }

  return(res)
}

# }}}

# extractPerformance {{{

#' @describeIn writePerformance Subset a performance table to the rows
#'   belonging to management procedures matching a pattern, together with their
#'   corresponding OM baseline rows (\code{mp == ""}).
#' @keywords manip

extractPerformance <- function(dat, mp) {

  # TODO: PARSE multiple MPs and match each OM

  # ASSIGN to avoid column match 
  smp <- mp

  # FIND mps & om
  if(length(smp) > 1)
    sub <- dat[mp %in% smp]
  else
    sub <- dat[mp %like% smp]
  mps <- sub[, as.character(unique(mp))]
  oms <- sub[, as.character(unique(om))]

  # RETURN subset, om + mp
  return(dat[om %in% oms & mp %in% c("", mps)])
}
# }}}

# getOMPerformance {{{

#' @describeIn writePerformance Batch-read serialized OM files, trim each to
#'   \code{fy} with [FLCore::window()], compute [mse::performance()], and
#'   row-bind the results into a single table.
#' @keywords file

getOMPerformance <- function(path, pattern="*.rds", fy, ...) {
  return(rbindlist(lapply(list.files(path, pattern, full.names=TRUE), function(i)
    suppressWarnings(performance(window(readRDS(i)$om, end=fy), ...)[,
      om:=sub('.rds', '', basename(i))])[])))
}
# }}}

# getMSEPerformance {{{

#' @describeIn writePerformance Batch-read serialized MSE result files or
#'   pre-computed performance \code{data.table} files and row-bind them into a
#'   single table.
#' @keywords file

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

#' Convert a performance statistics table to FLQuant objects
#'
#' Transforms a performance statistics `data.table`, or an `FLmse`/`FLmses`
#' object, into one or more [FLCore::FLQuant-class] objects organized by
#' management procedure, with statistics in the first (`quant`) dimension and
#' iterations preserved in the `iter` dimension.
#'
#' `performanceFLQuants()` always returns a named list (one `FLQuant` per MP).
#' `performanceFLQuant()` is a convenience wrapper that unwraps the list when
#' only one MP is present.
#'
#' @param x An object of class [mse::FLmse-class], [mse::FLmses-class], or a
#'   `data.table`. When `x` is an `FLmse` or `FLmses`, the performance slot is
#'   extracted automatically via [mse::performance()]. When `x` is a
#'   `data.table`, it must contain at least the columns `statistic`, `year`,
#'   `iter`, `data`, and `mp`. An optional `biol` column, present in outputs
#'   from [mse::mp()] on `FLombf` operating models, is mapped to the `unit`
#'   dimension of the resulting [FLCore::FLQuant-class].
#'
#' @return
#' \describe{
#'   \item{`performanceFLQuants()`}{A named list of [FLCore::FLQuant-class]
#'     objects, one per unique `mp` value. The first (`quant`) dimension is
#'     named by `statistic`; `year` and `iter` dimensions correspond to
#'     simulation years and iterations. When a `biol` column is present its
#'     values populate the `unit` dimension.}
#'   \item{`performanceFLQuant()`}{If the list contains exactly one element,
#'     that [FLCore::FLQuant-class] is returned directly; otherwise the full
#'     named list is returned.}
#' }
#'
#' @details
#' Conversion from `data.table` to `FLQuant` is delegated to
#' [FLCore::as.FLQuant()], which expects the columns to map unambiguously onto
#' the standard FLQuant dimensions. The `statistic` column populates the
#' `quant` dimension and must therefore contain values meaningful as dimension
#' names.
#'
#' When a `biol` column is present, [data.table::setnames()] renames it to
#' `unit` before passing to [FLCore::as.FLQuant()], leaving the original table
#' unmodified.
#'
#' An informative error is raised if `x` is a `data.table` that does not
#' contain all required columns (`statistic`, `year`, `iter`, `data`, `mp`).
#'
#' @seealso [mse::performance()], [mse::FLmse-class], [mse::FLmses-class],
#'   [FLCore::FLQuant-class], [FLCore::as.FLQuant()], [readPerformance()],
#'   [writePerformance()]
#'
#' @examples
#' \dontrun{
#' # From an FLmse object
#' flqs <- performanceFLQuants(mse_run)
#'
#' # From a pre-read data.table (multiple MPs)
#' dat <- readPerformance()
#' flqs <- performanceFLQuants(dat)
#' flqs[["MP1"]]["SSB", , , , ,]
#'
#' # Convenience wrapper, returns FLQuant directly for a single MP
#' flq <- performanceFLQuant(dat)
#' }
#'
#' @author Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#' @name performanceFLQuants
#' @keywords manip

performanceFLQuant <- function(x) {
  
  res <- performanceFLQuants(x)

  if(is.list(res) && length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}

#' @describeIn performanceFLQuants Convenience wrapper returning a single
#'   [FLCore::FLQuant-class] directly when only one management procedure is
#'   present, or the full named list otherwise.

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
