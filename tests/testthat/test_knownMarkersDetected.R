context("knownMarkersDetected")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seuratAllMarkers.rda"),
    package = "bcbioSingleCell"))
knownMarkers <- cellTypeMarkers[["hsapiens"]]

test_that("knownMarkersDetected", {
    data <- knownMarkersDetected(
        all = seuratAllMarkers,
        known = knownMarkers
    )
    expect_is(data, "grouped_df")
    group <- group_vars(data)
    expect_identical(group, "cell")

    annotable <- annotable(seurat)
    colnames(annotable)

    # Need better way to test P values here than using round.
    # Check old code, I think there's a stringr method that works well.
    # Or we can using `format()`.
    subset <- data[1L, setdiff(colnames(data), colnames(annotable))] %>%
        mutate_if(is.numeric, funs(round(., digits = 3L)))
    tbl <- tibble(
        "cell" = "Natural Killer Cell",
        "cluster" = factor("1", levels = c("0", "1", "2", "3")),
        pct1 = 0.613,
        pct2 = 0.862,
        avgLogFC = -1.018,
        pvalue = 0,
        padj = 0
    )
})
