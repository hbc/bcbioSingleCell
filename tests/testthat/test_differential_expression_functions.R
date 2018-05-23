context("Differential Expression Functions")



# diffExp ======================================================================
object <- cellranger_small
numerator <- colnames(object)[which(object[["sampleName"]] == "proximal")]
denominator <- colnames(object)[which(object[["sampleName"]] == "distal")]

test_that("diffExp : zinbwave-DESeq2", {
    x <- diffExp(
        object = object,
        numerator = numerator,
        denominator = denominator,
        zeroWeights = "zinbwave",
        caller = "DESeq2"
    )
    expect_s4_class(x, "DESeqResults")
})

test_that("diffExp : zinbwave-edgeR", {
    x <- diffExp(
        object = object,
        numerator = numerator,
        denominator = denominator,
        zeroWeights = "zinbwave",
        caller = "edgeR"
    )
    expect_s4_class(x, "DGELRT")
})

test_that("diffExp : zingeR-DESeq2", {
    x <- diffExp(
        object = object,
        numerator = numerator,
        denominator = denominator,
        zeroWeights = "zingeR",
        caller = "DESeq2"
    )
    expect_s4_class(x, "DESeqResults")
})

test_that("diffExp : zinbwave-edgeR", {
    x <- diffExp(
        object = object,
        numerator = numerator,
        denominator = denominator,
        zeroWeights = "zinbwave",
        caller = "edgeR"
    )
    expect_s4_class(x, "DGELRT")
})

test_that("diffExp : seurat", {
    # Expression in cluster 3 relative to cluster 2
    object <- Seurat::pbmc_small
    numerator <- Seurat::WhichCells(object, ident = 3L)
    denominator <- Seurat::WhichCells(object, ident = 2L)
    x <- diffExp(
        object = object,
        numerator = numerator,
        denominator = denominator,
        minCellsPerGene = 5L,
        minCountsPerCell = 5L
    )
    expect_s4_class(x, "DGELRT")
})
