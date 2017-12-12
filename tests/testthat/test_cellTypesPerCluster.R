context("cellTypesPerCluster")

load(system.file(
    file.path("extdata", "knownMarkersDetected.rda"),
    package = "bcbioSingleCell"))

test_that("cellTypesPerCluster", {
    data <- cellTypesPerCluster(knownMarkersDetected)
    expect_is(data, "grouped_df")

    group <- group_vars(data)
    expect_identical(group, "cluster")

    tbl <- tibble(
        cluster = factor("2", levels = c("0", "1", "2", "3")),
        cell = "Natural Killer Cell",
        n = 1L,
        symbol = "NCAM1",
        ensgene = "ENSG00000149294"
    ) %>%
        group_by(!!sym(group))
    expect_identical(data, tbl)
})
