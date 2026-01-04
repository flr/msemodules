# performance.R - DESC
# /home/mosqu003/Projects/FLR/code/mse/msemodules/msemodules/R/performance.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


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

#' Extracts performance time series for an MP including the corresponding historical OM
#'

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

# getPerformance

getOMPerformance <- function(path, pattern="*.rds", fy, ...) {
  return(rbindlist(lapply(list.files(path, pattern, full.names=TRUE), function(i)
    suppressWarnings(performance(window(readRDS(i)$om, end=fy), ...)[,
      om:=sub('.rds', '', basename(i))])[])))
}

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
