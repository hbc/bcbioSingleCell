set.seed(42L)

library(testthat)
library(bcbioSingleCell)

loadRemoteData(
    file.path(
        "http://bcbiosinglecell.seq.cloud",
        "testthat",
        "bcb.rda"),
    quiet = TRUE)

pooled <- suppressMessages(aggregateReplicates(bcb))
filtered <- filterCells(pooled, quiet = TRUE)

test_check("bcbioSingleCell")
