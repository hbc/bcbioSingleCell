context("prepareSingleCellTemplate")

test_that("prepareSingleCellTemplate", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "bibliography.bib",
        "setup.R")
    expect_silent(prepareSingleCellTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})
