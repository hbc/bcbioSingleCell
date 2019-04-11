context("plotUMIsPerCell")

point <- eval(formals(plotUMIsPerCell.bcbioSingleCell)[["point"]])
with_parameters_test_that(
    "ECDF inflection, knee points", {
        x <- plotUMIsPerCell(
            object = indrops,
            geom = "ecdf",
            point = point
        )
        expect_s3_class(x, "ggplot")
    },
    point = point
)
