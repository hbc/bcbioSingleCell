context("Plot Functions")



# plotCellCounts ===============================================================
test_that("plotCellCounts : bcbioSingleCell", {
    p <- plotCellCounts(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotCellCounts : seurat", {
    p <- plotCellCounts(seurat_small)
    expect_is(p, "ggplot")
})



# plotCellTypesPerCluster ======================================================
test_that("plotCellTypesPerCluster : seurat", {
    x <- cellTypesPerCluster(known_markers_small)[1L, , drop = FALSE]
    expect_identical(x[["geneName"]], "NCAM1")
    invisible(capture.output(
        p <- plotCellTypesPerCluster(
            object = seurat_small,
            cellTypesPerCluster = x
        )
    ))
    expect_is(p, "list")
    expect_is(p[[1L]][[1L]], "ggplot")
})



# plotGenesPerCell =============================================================
test_that("plotGenesPerCell : bcbioSingleCell", {
    p <- plotGenesPerCell(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotGenesPerCell : seurat", {
    p <- plotGenesPerCell(seurat_small)
    expect_is(p, "ggplot")
})



# plotFeatureTSNE ==============================================================
test_that("plotFeatureTSNE : seurat", {
    p <- plotFeatureTSNE(seurat_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})



# plotKnownMarkersDetected =====================================================
test_that("plotKnownMarkersDetected : seurat", {
    invisible(capture.output(
        p <- plotKnownMarkersDetected(
            object = seurat_small,
            markers = known_markers_small
        )
    ))
    expect_is(p, "list")
})



# plotMarkerTSNE ===============================================================
test_that("plotMarkerTSNE : seurat", {
    genes <- rownames(counts(Seurat::pbmc_small))[[1L]]
    args <- methodFormals("plotMarkerTSNE", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(args, function(arg) {
        p <- plotMarkerTSNE(Seurat::pbmc_small, genes = genes, expression = arg)
        expect_is(p, "ggplot")
    }))
})



# plotMarkers ==================================================================
test_that("plotMarkers : seurat", {
    genes <- head(rownames(seurat_small), n = 2L)
    invisible(capture.output(
        plotlist <- plotMarkers(
            object = seurat_small,
            genes = genes
        )
    ))
    expect_is(plotlist, "list")
    expect_identical(length(plotlist), 2L)
    expect_identical(names(plotlist), names(genes))
    p <- plotlist[[1L]]
    expect_is(p, "ggplot")
})



# plotMitoRatio ================================================================
test_that("plotMitoRatio : bcbioSingleCell", {
    p <- plotMitoRatio(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotMitoRatio : seurat", {
    p <- plotMitoRatio(seurat_small)
    expect_is(p, "ggplot")
})



# plotMitoVsCoding =============================================================
test_that("plotMitoVsCoding : bcbioSingleCell", {
    p <- plotMitoVsCoding(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotMitoVsCoding : seurat", {
    p <- plotMitoVsCoding(seurat_small)
    expect_is(p, "ggplot")
})



# plotNovelty ==================================================================
test_that("plotNovelty : bcbioSingleCell", {
    p <- plotNovelty(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotNovelty : seurat", {
    p <- plotNovelty(seurat_small)
    expect_is(p, "ggplot")
})



# plotPCA ======================================================================
test_that("plotPCA : seurat", {
    p <- plotPCA(seurat_small)
    expect_is(p, "ggplot")
})



# plotPCElbow ==================================================================
test_that("plotPCElbow : seurat", {
    x <- plotPCElbow(seurat_small)
    expect_identical(x, seq_len(9L))

    x <- plotPCElbow(Seurat::pbmc_small)
    expect_identical(x, seq_len(11L))
})



# plotQC =======================================================================
test_that("plotQC : grid", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    p <- plotQC(bcb_small, return = "grid")
    expect_is(p, "ggplot")
})

test_that("plotQC : list", {
    p <- plotQC(bcb_small, return = "list")
    expect_identical(
        names(p),
        c("plotReadsPerCell",
          "plotCellCounts",
          "plotUMIsPerCell",
          "plotGenesPerCell",
          "plotUMIsVsGenes",
          "plotMitoRatio",
          "plotNovelty")
    )
})

test_that("plotQC : markdown", {
    output <- capture.output(plotQC(bcb_small, return = "markdown"))
    sep <- c("", "", "")
    expect_identical(
        head(output, 3L),
        c("", "", "## Filtered quality control metrics {.tabset}")
    )
})

test_that("plotQC : seurat", {
    p <- plotQC(seurat_small)
    expect_is(p, "ggplot")
})



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
        plotReadsPerCell(Seurat::pbmc_small),
        "object does not contain nCount column in `metrics\\(\\)`"
    )
    expect_is(
        suppressWarnings(plotReadsPerCell(Seurat::pbmc_small)),
        "NULL"
    )
})



# plotTopMarkers ===============================================================
test_that("plotTopMarkers : seurat", {
    top <- topMarkers(all_markers_small, n = 1L)
    expect_is(top, "grouped_df")
    expect_identical(nrow(top), 5L)
    top <- top[seq_len(2L), ]
    invisible(capture.output(
        x <- plotTopMarkers(seurat_small, topMarkers = top)
    ))
    expect_is(x, "list")
    expect_is(x[[1L]][[1L]], "ggplot")
})



# plotTSNE =====================================================================
test_that("plotTSNE : seurat", {
    p <- plotTSNE(Seurat::pbmc_small)
    expect_is(p, "ggplot")
})



# plotUMIsPerCell ==============================================================
test_that("plotUMIsPerCell : bcbioSingleCell", {
    p <- plotUMIsPerCell(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotUMIsPerCell : seurat", {
    p <- plotUMIsPerCell(seurat_small)
    expect_is(p, "ggplot")
})



# plotUMIsVsGenes ==============================================================
test_that("plotUMIsVsGenes : bcbioSingleCell", {
    p <- plotUMIsVsGenes(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotUMIsVsGenes : seurat", {
    p <- plotUMIsVsGenes(seurat_small)
    expect_is(p, "ggplot")

    p <- plotUMIsVsGenes(Seurat::pbmc_small)
    expect_is(p, "ggplot")
})



# plotZerosVsDepth =============================================================
test_that("plotZerosVsDepth : bcbioSingleCell", {
    p <- plotZerosVsDepth(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotZerosVsDepth : seurat", {
    p <- plotZerosVsDepth(Seurat::pbmc_small)
    expect_is(p, "ggplot")
})
