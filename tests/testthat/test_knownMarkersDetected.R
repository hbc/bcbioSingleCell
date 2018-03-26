context("knownMarkersDetected")

test_that("knownMarkersDetected", {
    x <- knownMarkersDetected(
        all = all_markers_small,
        known = cellTypeMarkers[["homoSapiens"]]
    )
    expect_is(x, "grouped_df")
    group <- dplyr::group_vars(x)
    expect_identical(group, "cellType")

    annotable <- annotable(seurat_small)
    colnames(annotable)

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
