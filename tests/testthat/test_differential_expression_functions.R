context("Differential Expression Functions")



# diffExp ======================================================================
# Expression in cluster 3 relative to cluster 2
object <- seurat_small
numerator <- Seurat::WhichCells(object, ident = 3L)
denominator <- Seurat::WhichCells(object, ident = 2L)

test_that("diffExp : zinbwave-edgeR", {
    x <- diffExp(
        object = object,
        numerator = numerator,
        denominator = denominator,
        caller = "edgeR"
    )
    expect_s4_class(x, "DGELRT")
})

# DESeq2 is still relatively slow
test_that("diffExp : zinbwave-DESeq2", {
    x <- diffExp(
        object = object,
        numerator = numerator,
        denominator = denominator,
        caller = "DESeq2"
    )
    expect_s4_class(x, "DESeqResults")
})
