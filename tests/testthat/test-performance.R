# test-performance-db.R - DESC
# /home/mosqu003/Active/mse_FLR/msemodules/tests/testthat/test-performance-db.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


load('mse.rda')

# --- performanceFLQuant() {{{

tab <- performance(te)
res <- performanceFLQuant(te)

# CLASS
expect_true(is(res, "FLQuant"))

# DIM
tab[, .(statistic, year, iter, mp)]
dtab <- unname(unlist(tab[, lapply(.SD, uniqueN)])[c("statistic", "year", "iter")])

expect_equal(dim(res)[c(1,2,6)], dtab)

# }}}

# --- performance database helpers {{{

tmp <- tempfile(fileext=".dat.gz")
on.exit(unlink(tmp), add=TRUE)

dat <- data.table::copy(tab)
writePerformance(dat, file=tmp, overwrite=TRUE)

# writePerformance should not mutate caller object by reference
expect_true(is.numeric(dat$year))
expect_true(is.numeric(dat$data))

# writePerformance should use the mse performance schema
raw <- data.table::fread(tmp)
expect_identical(names(raw)[seq_along(mse:::.standardizeDT(data.table::copy(tab)))],
  names(mse:::.standardizeDT(data.table::copy(tab))))
expect_true(is.integer(raw$year))
expect_true(is.integer(raw$iter))
expect_true(is.numeric(raw$data))

# hasPerformance / listPerformance basic behaviour
expect_true(hasPerformance(file=tmp))
expect_true(nrow(listPerformance(file=tmp)) > 0)

# validatePerformance should fail before write-time coercion
bad <- data.table::copy(tab)
bad[, year := as.character(year)]
expect_error(writePerformance(bad, file=tmp, overwrite=TRUE),
  "Column 'year' must be numeric.")

# readPerformance should tolerate missing optional label column
tmp_nolab <- tempfile(fileext=".dat.gz")
on.exit(unlink(tmp_nolab), add=TRUE)
no_label <- data.table::copy(tab)[, label := NULL]
data.table::fwrite(no_label, file=tmp_nolab)
re <- readPerformance(tmp_nolab)
expect_false("label" %in% names(re))

# extractPerformance should support exact matching for vector mp
ex <- data.table::data.table(
  om=c("om1", "om1", "om1", "om2", "om2"),
  mp=c("", "hcrA", "hcrB", "", "hcrB")
)
sub <- extractPerformance(ex, c("hcrA", "hcrB"))
expect_equal(sort(unique(as.character(sub$mp))), c("", "hcrA", "hcrB"))
expect_equal(sort(unique(as.character(sub$om))), c("om1", "om2"))

# diffPerformance and deletePerformance
mod <- data.table::copy(dat)
mod[1, data := data + 1]
d <- diffPerformance(mod, file=tmp)
expect_equal(sort(names(d)), c("new", "replace", "unchanged"))
expect_true(data.table::is.data.table(d$new))
expect_true(data.table::is.data.table(d$replace))
expect_true(data.table::is.data.table(d$unchanged))

deleted <- deletePerformance(file=tmp, om=as.character(dat$om[1]), dry_run=TRUE)
expect_true(nrow(deleted) > 0)
expect_true(file.exists(tmp))

# }}}
