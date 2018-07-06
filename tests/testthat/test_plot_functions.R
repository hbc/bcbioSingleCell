context("Plot Functions")



# plotCellCounts ===============================================================
test_that("plotCellCounts : bcbioSingleCell", {
    p <- plotCellCounts(indrops_small)
    expect_is(p, "ggplot")
})

test_that("plotCellCounts : seurat", {
    p <- plotCellCounts(seurat_small)
    expect_is(p, "ggplot")
})



# plotCellTypesPerCluster ======================================================
test_that("plotCellTypesPerCluster : seurat", {
    cellTypesPerCluster <- cellTypesPerCluster(known_markers_small) %>%
        # Subset for speed
        head(2L)
    invisible(capture.output(
        p <- plotCellTypesPerCluster(
            object = seurat_small,
            cellTypesPerCluster = cellTypesPerCluster
        )
    ))
    expect_is(p, "list")
    expect_is(p[[1L]][[1L]], "ggplot")
})



# plotGenesPerCell =============================================================
test_that("plotGenesPerCell : bcbioSingleCell", {
    p <- plotGenesPerCell(indrops_small)
    expect_is(p, "ggplot")
})

test_that("plotGenesPerCell : seurat", {
    p <- plotGenesPerCell(seurat_small)
    expect_is(p, "ggplot")
})



# plotFeature ==================================================================
test_that("plotFeatureTSNE : seurat", {
    p <- plotFeatureTSNE(seurat_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})

test_that("plotFeatureUMAP : seurat", {
    p <- plotFeatureUMAP(seurat_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})



# plotMarker ===================================================================
object <- seurat_small
genes <- head(rownames(object))

test_that("plotMarkerTSNE : seurat", {
    args <- methodFormals("plotMarkerTSNE", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(args, function(arg) {
        p <- plotMarkerTSNE(object, genes = genes, expression = arg)
        expect_is(p, "ggplot")
    }))
})

test_that("plotMarkerUMAP : seurat", {
    args <- methodFormals("plotMarkerUMAP", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(args, function(arg) {
        p <- plotMarkerUMAP(object, genes = genes, expression = arg)
        expect_is(p, "ggplot")
    }))
})

test_that("plotKnownMarkersDetected : seurat", {
    invisible(capture.output(
        p <- plotKnownMarkersDetected(
            object = seurat_small,
            markers = head(known_markers_small, 2L)
        )
    ))
    expect_is(p, "list")
})

test_that("plotTopMarkers : seurat", {
    markers <- topMarkers(all_markers_small, n = 1L) %>%
        # Subset for speed
        head(2L)
    invisible(capture.output(
        x <- plotTopMarkers(
            object = seurat_small,
            markers = markers
        )
    ))
    expect_is(x, "list")
    expect_is(x[[1L]][[1L]], "ggplot")
})



# plotMitoRatio ================================================================
test_that("plotMitoRatio : bcbioSingleCell", {
    p <- plotMitoRatio(indrops_small)
    expect_is(p, "ggplot")
})

test_that("plotMitoRatio : seurat", {
    p <- plotMitoRatio(seurat_small)
    expect_is(p, "ggplot")
})



# plotMitoVsCoding =============================================================
test_that("plotMitoVsCoding : bcbioSingleCell", {
    p <- plotMitoVsCoding(indrops_small)
    expect_is(p, "ggplot")
})

test_that("plotMitoVsCoding : seurat", {
    p <- plotMitoVsCoding(seurat_small)
    expect_is(p, "ggplot")
})



# plotNovelty ==================================================================
test_that("plotNovelty : bcbioSingleCell", {
    p <- plotNovelty(indrops_small)
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
    expect_identical(x, seq_len(7L))

    x <- plotPCElbow(Seurat::pbmc_small)
    expect_identical(x, seq_len(11L))
})



# plotQC =======================================================================
test_that("plotQC : grid", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    p <- plotQC(indrops_small, return = "grid")
    expect_is(p, "ggplot")
})

test_that("plotQC : list", {
    p <- plotQC(indrops_small, return = "list")
    expect_identical(
        names(p),
        c(
            "Cell Counts",
            "Reads per Cell",
            "UMIs per Cell",
            "Genes per Cell",
            "UMIs vs. Genes",
            "Novelty",
            "Mito Ratio",
            "Zeros vs. Depth"
        )
    )
})

test_that("plotQC : markdown", {
    output <- capture.output(plotQC(indrops_small, return = "markdown"))
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
# Example dataset doesn't have a cellular barcode cutoff because we removed the
# bcbio commands log file (which conflicts with Travis CI)
test_that("plotReadsPerCell : bcbioSingleCell", {
    # Histogram
    x <- plotReadsPerCell(indrops_small, geom = "histogram")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomStep"
    )

    # Ridgeline
    x <- plotReadsPerCell(indrops_small, geom = "ridgeline")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomDensityRidges"
    )

    # Violin
    x <- plotReadsPerCell(indrops_small, geom = "violin")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomViolin"
    )

    # ECDF
    x <- plotReadsPerCell(indrops_small, geom = "ecdf")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomStep"
    )
})



# plotTSNE =====================================================================
test_that("plotTSNE : seurat", {
    p <- plotTSNE(seurat_small)
    expect_is(p, "ggplot")
})



# plotUMAP =====================================================================
test_that("plotUMAP : seurat", {
    p <- plotTSNE(seurat_small)
    expect_is(p, "ggplot")
})



# plotUMIsPerCell ==============================================================
test_that("plotUMIsPerCell : bcbioSingleCell", {
    p <- plotUMIsPerCell(indrops_small)
    expect_is(p, "ggplot")
})

test_that("plotUMIsPerCell : seurat", {
    p <- plotUMIsPerCell(seurat_small)
    expect_is(p, "ggplot")
})



# plotUMIsVsGenes ==============================================================
test_that("plotUMIsVsGenes : bcbioSingleCell", {
    p <- plotUMIsVsGenes(indrops_small)
    expect_is(p, "ggplot")
})

test_that("plotUMIsVsGenes : seurat", {
    p <- plotUMIsVsGenes(seurat_small)
    expect_is(p, "ggplot")

    p <- plotUMIsVsGenes(seurat_small)
    expect_is(p, "ggplot")
})



# plotZerosVsDepth =============================================================
test_that("plotZerosVsDepth : bcbioSingleCell", {
    p <- plotZerosVsDepth(indrops_small)
    expect_is(p, "ggplot")
})

test_that("plotZerosVsDepth : seurat", {
    p <- plotZerosVsDepth(seurat_small)
    expect_is(p, "ggplot")
})
