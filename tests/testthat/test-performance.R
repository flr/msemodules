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

