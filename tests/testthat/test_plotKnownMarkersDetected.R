context("plotKnownMarkersDetected")

test_that("plotFeatureTSNE", {
    p <- plotFeatureTSNE(seurat, features = c("PC1", "PC2"))
    expect_is(p, "ggplot")
})
