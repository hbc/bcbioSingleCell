context("plotUMIsPerCell")

skip_if_not(packageVersion("DropletUtils") >= "1.4")
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
