context("plotReadsPerCell")

# Example dataset doesn't have a cellular barcode cutoff because we removed the
# bcbio commands log file (which conflicts with Travis CI).

geom <- eval(formals(plotReadsPerCell.bcbioSingleCell)[["geom"]])
with_parameters_test_that(
    "geom", {
        x <- plotReadsPerCell(indrops, geom = geom)
        expect_is(x, "ggplot")
    },
    geom = geom
)
