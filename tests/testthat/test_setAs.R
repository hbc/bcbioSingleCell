context("setAs coercion")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "filtered.rda"),
    package = "bcbioSingleCell"))

test_that("Coerce bcbioSingleCell to seurat", {
    seurat <- as(filtered, "seurat")
    expect_is(seurat, "seurat")

    # Check slotted count integrity
    counts <- slot(seurat, "raw.data")
    expect_is(counts, "dgCMatrix")
    expect_identical(
        dim(counts),
        c(1000L, 243L)
    )

    # Check to make sure our Seurat coercion step isn't dropping data
    origcounts <- assay(filtered)
    expect_identical(
        dim(counts),
        dim(origcounts)
    )
})

test_that("Expected warnings", {
    expect_warning(
        as(bcb, "seurat"),
        "Use `aggregateReplicates\\(\\)` to merge technical replicates"
    )
    expect_warning(
        as(bcb, "seurat"),
        "Use `filterCells\\(\\)` to apply filtering cutoffs"
    )
})
