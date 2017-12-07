set.seed(42L)

library(testthat)
library(bcbioSingleCell)

loadRemoteData(
    file.path(
        "http://bcbiosinglecell.seq.cloud",
        "testthat",
        "bcb.rda"),
    quiet = TRUE)

test_check("bcbioSingleCell")
