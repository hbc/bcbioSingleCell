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
        object = all_markers_small,
        known = cell_type_markers[["homoSapiens"]]
    )
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "cellType")
    expect_identical(
        lapply(x, class),
        list(
            cellType = "factor",
            cluster = "factor",
            geneID = "character",
            geneName = "character",
            avgLogFC = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            rowname = "character",
            pvalue = "numeric",
            padj = "numeric",
            geneBiotype = "factor",
            description = "character",
            seqCoordSystem = "factor",
            broadClass = "factor"
        )
    )
})



# sanitizeMarkers ==============================================================
test_that("sanitizeMarkers : seurat", {
    expect_message(
        sanitizeMarkers(
            object = seurat_small,
            markers = all_markers_small
        ),
        "Markers are already sanitized"
    )

    # FindAllMarkers
    invisible(capture.output(
        all <- Seurat::FindAllMarkers(seurat_small)
    ))
    x <- sanitizeMarkers(object = seurat_small, markers = all)
    expect_is(x, "grouped_df")

    # FindMarkers
    invisible(capture.output(
        ident3 <- Seurat::FindMarkers(
            seurat_small,
            ident.1 = "3",
            ident.2 = NULL
        )
    ))
    x <- sanitizeMarkers(
        object = seurat_small,
        markers = ident3
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
        lapply(x, class),
        list(
            cluster = "factor",
            avgLogFC = "numeric",
            pct1 = "numeric",
            pct2 = "numeric",
            rowname = "character",
            pvalue = "numeric",
            padj = "numeric",
            geneID = "character",
            geneName = "character",
            geneBiotype = "factor",
            description = "character",
            seqCoordSystem = "factor",
            broadClass = "factor"
        )
    )
})
