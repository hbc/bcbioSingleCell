context("Clustering Functions")



# cellTypesPerCluster ==========================================================
test_that("cellTypesPerCluster", {
    x <- cellTypesPerCluster(known_markers_small)
    expect_is(x, "grouped_df")
    group <- dplyr::group_vars(x)
    expect_identical(group, "cluster")
    y <- tibble::tibble(
        "cluster" = factor("1", levels = c("0", "1", "2", "3", "4")),
        "cellType" = "Natural Killer Cell",
        "n" = 1L,
        "geneID" = "ENSG00000149294",
        "geneName" = "NCAM1"
    ) %>%
        dplyr::group_by(!!rlang::sym(group))
    expect_identical(x, y)
})



# knownMarkersDetected =========================================================
test_that("knownMarkersDetected", {
    x <- knownMarkersDetected(
        object = all_markers_small,
        known = cellTypeMarkers[["homoSapiens"]]
    )
    expect_is(x, "grouped_df")
    group <- dplyr::group_vars(x)
    expect_identical(group, "cellType")

    # Need better way to test P values here than using round.
    # Check old code, I think there's a stringr method that works well.
    # Or we can using `format()`.
    subset <- x %>%
        .[1L, , drop = FALSE] %>%
        dplyr::mutate_if(
            is.numeric,
            dplyr::funs(round(., digits = 3L))
        )
    seqnamesLevels <- levels(subset[["seqnames"]])
    strandLevels <- levels(subset[["strand"]])
    geneBiotypeLevels <- levels(subset[["geneBiotype"]])
    seqCoordSystemLevels <- levels(subset[["seqCoordSystem"]])
    broadClassLevels <- levels(subset[["broadClass"]])
    target <- tibble::tibble(
        "cellType" = "Natural Killer Cell",
        "cluster" = factor("1", levels = c("0", "1", "2", "3", "4")),
        "geneID" = "ENSG00000149294",
        "geneName" = "NCAM1",
        "pct1" = 0.466,
        "pct2" = 0.701,
        "avgLogFC" = -0.954,
        "pvalue" = 0,
        "padj" = 0,
        "rowname" = "NCAM1",
        "seqnames" = factor("11", levels = seqnamesLevels),
        "start" = 112961247,
        "end" = 113278436,
        "width" = 317190,
        "strand" = factor("+", levels = strandLevels),
        "geneBiotype" = factor("protein_coding", levels = geneBiotypeLevels),
        "description" = paste(
            "neural cell adhesion molecule 1",
            "[Source:HGNC Symbol;Acc:HGNC:7656]"
        ),
        "seqCoordSystem" = factor("chromosome", levels = seqCoordSystemLevels),
        "broadClass" = factor("coding", levels = broadClassLevels)
    )
    expect_equal(subset, target)
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
    all <- FindAllMarkers(seurat_small)
    x <- sanitizeMarkers(
        object = seurat_small,
        markers = all
    )
    expect_is(x, "grouped_df")

    # FindMarkers
    ident3 <- FindMarkers(seurat_small, ident.1 = "3", ident.2 = NULL)
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
    expect_identical(
        dplyr::group_vars(x),
        "cluster"
    )
    expect_identical(
        dim(x),
        c(50L, 18L)
    )
})
