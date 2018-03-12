context("cellTypesPerCluster")

test_that("cellTypesPerCluster", {
    x <- cellTypesPerCluster(known_markers_detected)
    expect_is(x, "grouped_df")
    group <- group_vars(x)
    expect_identical(group, "cluster")
    # TODO Switch to ident here from cluster
    # TODO Switch to `cellType`
    tbl <- tibble(
        "cluster" = factor("2", levels = c("0", "1", "2", "3")),
        "cell" = "Natural Killer Cell",
        "n" = 1L,
        "geneID" = "ENSG00000149294",
        "geneName" = "NCAM1"
    ) %>%
        group_by(!!sym(group))
    expect_identical(x, tbl)
})
