context("Clustering Functions")



# cellTypesPerCluster ==========================================================
test_that("cellTypesPerCluster", {
    x <- cellTypesPerCluster(known_markers_small)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cluster")
    expect_identical(
        lapply(x, class),
        list(
            cluster = "factor",
            cellType = "factor",
            n = "integer",
            geneID = "character",
            geneName = "character"
        )
    )
})



# knownMarkersDetected =========================================================
test_that("knownMarkersDetected", {
    x <- knownMarkersDetected(
        all = all_markers_small,
        known = cell_type_markers[["homoSapiens"]]
    )
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cellType")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            avgLogFC = "numeric",
            broadClass = "factor",
            cellType = "factor",
            cluster = "factor",
            description = "factor",
            end = "integer",
            geneBiotype = "factor",
            geneID = "character",
            geneName = "factor",
            padj = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            pvalue = "numeric",
            rowname = "character",
            seqCoordSystem = "factor",
            seqnames = "factor",
            start = "integer",
            strand = "factor",
            width = "integer"
        )
    )
})



# sanitizeSeuratMarkers ========================================================
test_that("sanitizeSeuratMarkers", {
    # Early return on sanitized data
    expect_message(
        sanitizeSeuratMarkers(
            data = all_markers_small,
            rowRanges = rowRanges(seurat_small)
        ),
        "Markers are already sanitized"
    )

    # FindAllMarkers
    invisible(capture.output(
        all <- Seurat::FindAllMarkers(seurat_small)
    ))
    x <- sanitizeSeuratMarkers(
        data = all,
        rowRanges = rowRanges(seurat_small)
    )
    expect_is(x, "grouped_df")

    # FindMarkers
    invisible(capture.output(
        ident3 <- Seurat::FindMarkers(
            seurat_small,
            ident.1 = "3",
            ident.2 = NULL
        )
    ))
    x <- sanitizeSeuratMarkers(
        data = ident3,
        rowRanges = rowRanges(seurat_small)
    )
    expect_is(x, "data.frame")
    expect_true(tibble::has_rownames(x))
})



# topMarkers ===================================================================
test_that("topMarkers : grouped_df", {
    x <- topMarkers(all_markers_small)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cluster")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            avgLogFC = "numeric",
            broadClass = "factor",
            cluster = "factor",
            description = "factor",
            end = "integer",
            geneBiotype = "factor",
            geneID = "character",
            geneName = "factor",
            padj = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            pvalue = "numeric",
            rowname = "character",
            seqCoordSystem = "factor",
            seqnames = "factor",
            start = "integer",
            strand = "factor",
            width = "integer"
        )
    )
})
