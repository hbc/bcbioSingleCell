context("plotFeatureTSNE")

load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

test_that("plotFeatureTSNE", {
    p <- plotFeatureTSNE(seurat, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})
