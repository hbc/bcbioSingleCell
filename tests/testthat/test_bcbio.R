context("bcbio")

test_that("bcbio : seurat_small", {
    x <- bcbio(seurat_small)
    expect_is(x, "list")
    expect_identical(
        lapply(x, class),
        list(
            "rowRanges" = structure("GRanges", package = "GenomicRanges"),
            "metadata" = "list"
        )
    )
})

test_that("bcbio : seurat_small assignment", {
    x <- seurat_small
    # Stash as new slot
    bcbio(x, "stash") <- "XXX"
    expect_identical(
        bcbio(x, "stash"),
        "XXX"
    )
    # Metadata stash
    bcbio(x, "metadata")[["stash"]] <- "YYY"
    expect_identical(
        bcbio(x, "metadata")[["stash"]],
        "YYY"
    )
})
