context("plotPCElbow")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("seurat", {
    pcUse <- plotPCElbow(seurat)
    expect_identical(
        pcUse,
        seq_len(10L)
    )
})
