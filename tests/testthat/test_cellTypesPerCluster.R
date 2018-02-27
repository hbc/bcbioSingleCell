context("cellTypesPerCluster")

test_that("cellTypesPerCluster", {
    x <- cellTypesPerCluster(known_markers_detected)
    expect_is(x, "grouped_df")
    group <- group_vars(dx)
    expect_identical(group, "cluster")
    tbl <- tibble(
        cluster = factor("2", levels = c("0", "1", "2", "3")),
        cell = "Natural Killer Cell",
        n = 1L,
        symbol = "NCAM1",
        ensgene = "ENSG00000149294"
    ) %>%
        group_by(!!sym(group))
    expect_identical(x, tbl)
})
