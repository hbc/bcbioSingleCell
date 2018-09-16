context("Read Functions")



# readCellRanger ===============================================================
test_that("readCellRanger", {
    object <- suppressWarnings(readCellRanger(
        uploadDir = system.file(
            "extdata/cellranger",
            package = "bcbioSingleCell"
        ),
        organism = "Homo sapiens"
    ))
    expect_is(object, "CellRanger")
    expect_identical(dim(object), c(500L, 500L))
    expect_identical(
        sampleNames(object),
        c(pbmc_1 = "pbmc")
    )
    expect_identical(
        object = head(colnames(object), n = 1L),
        expected = "pbmc_1_AAACCTGAGAAGGCCT"
    )
    expect_identical(
        object = head(rownames(object), n = 1L),
        expected = "ENSG00000004487"
    )
})
