context("readCellTypeMarkers")

test_that("readCellTypeMarkers : Mus musculus", {
    file <- system.file(
        file.path("extdata", "cell_type_markers.csv"),
        package = "bcbioSingleCell"
    )
    gene2symbol <- makeGene2symbolFromEnsembl("Homo sapiens")
    x <- readCellTypeMarkers(
        file = file,
        gene2symbol = gene2symbol
    )
    group <- dplyr::group_vars(x)
    expect_identical(group, "cellType")
    x <- x[1L, ]
    y <- tibble::tibble(
        "cellType" = "B Cell",
        "geneID" = "ENSG00000177455",
        "geneName" = "CD19"
    ) %>%
        dplyr::group_by(!!rlang::sym("cellType"))
    expect_identical(x, y)
})
