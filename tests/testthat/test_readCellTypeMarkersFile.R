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
    expect_identical(group, "cell")
    expect_identical(
        data[1L, ],
        tibble(
            "cell" = "B Cell",
            "symbol" = "Blnk",
            "ensgene" = "ENSMUSG00000061132"
        ) %>%
            group_by(.data[["cell"]])
    )
})
