context("readCellTypeMarkers")

test_that("readCellTypeMarkers : Mus musculus", {
    file <- system.file(
        file.path("extdata", "cell_type_markers.csv"),
        package = "bcbioSingleCell"
    )
    gene2symbol <- makeGene2symbolFromEnsembl("Mus musculus")
    x <- readCellTypeMarkers(
        file = file,
        gene2symbol = gene2symbol
    )
    group <- dplyr::group_vars(x)
    expect_identical(group, "cellType")
    x <- x[1L, ]
    y <- tibble::tibble(
        "cellType" = "B Cell",
        "geneID" = "ENSMUSG00000061132",
        "geneName" = "Blnk"
    ) %>%
        dplyr::group_by(!!rlang::sym("cellType"))
    expect_identical(x, y)
})
