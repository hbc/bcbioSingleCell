context("gene2symbol")

colnames <- c("geneID", "geneName")

test_that("bcbioSingleCell", {
    x <- gene2symbol(bcb_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})

test_that("seurat_small", {
    x <- gene2symbol(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})
