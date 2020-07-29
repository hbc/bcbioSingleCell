context("plotReadsPerCell")

## Example dataset doesn't have a cellular barcode cutoff because we removed the
## bcbio commands log file (which conflicts with Travis CI).
test_that("geom", {
    for (geom in eval(formals(`plotReadsPerCell,bcbioSingleCell`)[["geom"]])) {
        x <- plotReadsPerCell(bcb, geom = geom)
        expect_is(x, "ggplot")
    }
})
