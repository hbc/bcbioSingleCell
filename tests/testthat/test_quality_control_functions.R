context("Quality Control Functions")



# plotReadsPerCell =============================================================
test_that("plotReadsPerCell : bcbioSingleCell", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    histogram <- plotReadsPerCell(bcb_small, geom = "histogram")
    expect_is(histogram, "ggplot")
    expect_is(
        histogram %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomLine"
    )

    ridgeline <- plotReadsPerCell(bcb_small, geom = "ridgeline")
    expect_is(ridgeline, "ggplot")
    expect_is(
        ridgeline %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomDensityRidges"
    )

    violin <- plotReadsPerCell(bcb_small, geom = "violin")
    expect_is(violin, "ggplot")
    expect_is(
        violin %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomViolin"
    )
})

test_that("plotReadsPerCell : seurat", {
    p <- plotReadsPerCell(seurat_small)
    expect_is(p, "ggplot")

    # seurat object not created by bcbioSingleCell
    expect_warning(
        plotReadsPerCell(pbmc_small),
        "object does not contain nCount column in `metrics\\(\\)`"
    )
})
