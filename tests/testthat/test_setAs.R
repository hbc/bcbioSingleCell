context("setAs coercion")

test_that("Coerce bcbioSingleCell to seurat", {
    seurat <- as(filtered, "seurat")
    expect_is(seurat, "seurat")

    # Check slotted count integrity
    counts <- slot(seurat, "raw.data")
    expect_is(counts, "dgCMatrix")
    expect_identical(
        dim(counts),
        c(21576L, 3292L)
    )

    # Check to make sure our Seurat coercion step isn't dropping data
    origcounts <- assay(filtered)
    expect_identical(
        dim(counts),
        dim(origcounts)
    )
})

test_that("Require technical replicate aggregation", {
    expect_error(
        as(bcb, "seurat"),
        paste("'aggregateReplicates\\(\\)' required to merge",
              "technical replicates prior to seurat coercion")
    )
})

test_that("Require filtered counts", {
    expect_error(
        as(pooled, "seurat"),
        "'filterCells\\(\\)' hasn't been applied to this dataset"
    )
})
