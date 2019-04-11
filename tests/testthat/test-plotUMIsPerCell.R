context("plotUMIsPerCell")

plotUMIsPerCell.bcbioSingleCell
point
with_parameters_test_that(
    "ECDF inflection, knee points", {
        x <- plotUMIsPerCell(
            object = indrops,
            geom = "ecdf",
            point = "inflection"
        )
        expect_s3_class(x, "ggplot")
    }

)
