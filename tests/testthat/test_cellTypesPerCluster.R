context("cellTypesPerCluster")

test_that("cellTypesPerCluster", {
    x <- cellTypesPerCluster(known_markers_small)
    expect_is(x, "grouped_df")
    group <- dplyr::group_vars(x)
    expect_identical(group, "cluster")
    y <- tibble::tibble(
        "cluster" = factor("1", levels = c("0", "1", "2", "3", "4")),
        "cellType" = "Natural Killer Cell",
        "n" = 1L,
        "geneID" = "ENSG00000149294",
        "geneName" = "NCAM1"
    ) %>%
        dplyr::group_by(!!sym(group))
    expect_identical(x, y)
})
