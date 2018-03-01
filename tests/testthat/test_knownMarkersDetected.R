context("knownMarkersDetected")

test_that("knownMarkersDetected", {
    known <- cellTypeMarkers[["homoSapiens"]]
    x <- knownMarkersDetected(
        all = seurat_all_markers,
        known = known)
    expect_is(x, "grouped_df")
    group <- group_vars(x)
    expect_identical(group, "cell")

    annotable <- annotable(seurat)
    colnames(annotable)

    # Need better way to test P values here than using round.
    # Check old code, I think there's a stringr method that works well.
    # Or we can using `format()`.
    subset <- x %>%
        .[1L, setdiff(colnames(.), colnames(annotable))] %>%
        mutate_if(is.numeric, funs(round(., digits = 3L)))
    tbl <- tibble(
        "cell" = "Natural Killer Cell",
        "cluster" = factor("1", levels = c("0", "1", "2", "3")),
        pct1 = 0.613,
        pct2 = 0.862,
        avgLogFC = -1.018,
        pvalue = 0L,
        padj = 0L
    )
})
