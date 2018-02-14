context("prepareSingleCellTemplate")

test_that("prepareSingleCellTemplate", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "_setup.R",
        "bibliography.bib"
    )
    expect_silent(prepareSingleCellTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})
