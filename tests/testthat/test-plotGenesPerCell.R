context("plotGenesPerCell")

geom <- eval(formals(`plotGenesPerCell,bcbioSingleCell`)[["geom"]])
with_parameters_test_that(
    "max/min argument", {
        p <- plotGenesPerCell(
            object = indrops,
            geom = geom,
            min = 1L,
            max = 100L
        )
        expect_s3_class(p, "ggplot")
    },
    geom = geom
)
