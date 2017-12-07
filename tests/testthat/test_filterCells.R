context("filterCells")

test_that("Default parameters", {
    filtered <- filterCells(bcb, quiet = TRUE)
    expect_is(filtered, "bcbioSingleCell")
    expect_is(
        metadata(filtered)[["filterParams"]],
        "numeric"
    )
    expect_is(
        metadata(filtered)[["filterCells"]],
        "character"
    )
    expect_is(
        metadata(filtered)[["filterGenes"]],
        "character"
    )
    expect_identical(
        metadata(filtered)[["subset"]],
        TRUE
    )
})

test_that("Capture output", {
    output <- capture.output(suppressMessages(
        filterCells(bcb, quiet = FALSE)
    ))
    expect_is(output, "character")
    expect_identical(
        output[[1]],
        "Filtering parameters:"
    )
})

test_that("Maximum parameters", {
    # This should return an object with the same dimensions
    filtered <- filterCells(
        bcb,
        minUMIs = 0,
        minGenes = 0,
        maxGenes = Inf,
        maxMitoRatio = 1,
        minNovelty = 0,
        minCellsPerGene = 0,
        quiet = TRUE)
    expect_is(filtered, "bcbioSingleCell")
    expect_equal(
        dim(filtered),
        dim(bcb)
    )
})
