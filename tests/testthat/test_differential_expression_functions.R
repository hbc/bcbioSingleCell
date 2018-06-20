context("Differential Expression Functions")



# differentialExpression =======================================================
# Expression in cluster 3 relative to cluster 2
object <- seurat_small
numerator <- Seurat::WhichCells(object, ident = 3L)
denominator <- Seurat::WhichCells(object, ident = 2L)

test_that("differentialExpression : zinbwave-edgeR", {
    x <- differentialExpression(
        object = object,
        numerator = numerator,
        denominator = denominator,
        zeroWeights = "zinbwave",
        caller = "edgeR"
    )
    expect_s4_class(x, "DGELRT")
})

test_that("differentialExpression : zingeR-edgeR", {
    x <- differentialExpression(
        object = object,
        numerator = numerator,
        denominator = denominator,
        zeroWeights = "zingeR",
        caller = "edgeR"
    )
    expect_s4_class(x, "DGELRT")
})

# DESeq2 is still relatively slow
test_that("differentialExpression : zinbwave-DESeq2", {
    x <- differentialExpression(
        object = object,
        numerator = numerator,
        denominator = denominator,
        zeroWeights = "zinbwave",
        caller = "DESeq2"
    )
    expect_s4_class(x, "DESeqResults")
})

# nolint start
# zingeR isn't importing DESeqDataSetFromMatrix correctly
# Filed an issue on GitHub, need to make a pull request
# test_that("differentialExpression : zingeR-DESeq2", {
#     x <- differentialExpression(
#         object = object,
#         numerator = numerator,
#         denominator = denominator,
#         zeroWeights = "zingeR",
#         caller = "DESeq2"
#     )
#     expect_s4_class(x, "DESeqResults")
# })
# nolint end
