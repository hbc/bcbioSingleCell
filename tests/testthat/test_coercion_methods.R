context("Coercion Methods")

test_that("Coerce bcbioSingleCell to seurat", {
    x <- as(indrops_small, "seurat")
    expect_is(x, "seurat")
    # Check slotted count integrity
    counts <- counts(x)
    expect_is(counts, "dgCMatrix")
    expect_identical(
        dim(counts),
        dim(indrops_small)
    )
})

# Require filtered counts
# Require technical replicate aggregation?
