context("plotKnownMarkersDetected")

test_that("plotFeatureTSNE", {
    p <- plotFeatureTSNE(seurat_small, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})
