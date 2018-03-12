context("readCellTypeMarkersFile")

test_that("Mus musculus", {
    cellTypeMarkersFile <- system.file(
        file.path("extdata", "cellTypeMarkers.csv"),
        package = "bcbioSingleCell")
    gene2symbol <- annotable("Mus musculus", format = "gene2symbol")
    data <- readCellTypeMarkersFile(
        cellTypeMarkersFile,
        gene2symbol = gene2symbol)
    group <- group_vars(data)
    expect_identical(group, "cellType")
    expect_identical(
        data[1L, ],
        # TODO Switch to cellType
        tibble(
            "cellType" = "B Cell",
            "geneID" = "ENSMUSG00000061132",
            "geneName" = "Blnk",
        ) %>%
            group_by(.data[["cellType"]])
    )
})
