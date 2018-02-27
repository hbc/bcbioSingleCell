context("plotReadsPerCell")

test_that("bcbioSingleCell", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    histogram <- plotReadsPerCell(bcb, geom = "histogram")
    expect_is(histogram, "ggplot")
    expect_is(
        histogram %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomLine"
    )

    ridgeline <- plotReadsPerCell(bcb, geom = "ridgeline")
    expect_is(ridgeline, "ggplot")
    expect_is(
        ridgeline %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomDensityRidges"
    )

    violin <- plotReadsPerCell(bcb, geom = "violin")
    expect_is(violin, "ggplot")
    expect_is(
        violin %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomViolin"
    )
})

test_that("filtered", {
    p <- plotReadsPerCell(filtered)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotReadsPerCell(seurat)
    expect_is(p, "NULL")
})
