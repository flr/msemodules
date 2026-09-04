# zzz.R - DESC
# /home/mosqu003/Projects/FLR/mse/msemodules/R/zzz.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

.onLoad <- function(libname, pkgname) {

  # register immediately if B is already installed
  if (requireNamespace("FLa4a", quietly = TRUE)) {
    .asFLoma4aFitSA()
  }

  # register lazily if B gets attached later
  setHook(
    packageEvent("FLa4a", "onLoad"),
    function(...) .asFLoma4aFitSA()
  )
}
