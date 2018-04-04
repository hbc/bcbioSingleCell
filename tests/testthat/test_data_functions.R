context("Data Functions")



# aggregateReplicates ==========================================================
test_that("aggregateReplicates", {
    x <- aggregateReplicates(bcb_small)
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(
        dim(x),
        c(500L, 500L)
    )
    map <- metadata(x)[["aggregateReplicates"]]
    expect_is(map, "factor")
    expect_identical(length(map), 500L)
    expect_identical(length(levels(map)), 500L)
})



# bcbio ========================================================================
test_that("bcbio : seurat_small", {
    x <- bcbio(seurat_small)
    expect_is(x, "list")
    expect_identical(
        lapply(x, class),
        list(
            "rowRanges" = structure("GRanges", package = "GenomicRanges"),
            "metadata" = "list"
        )
    )
})

test_that("bcbio : seurat_small assignment", {
    x <- seurat_small
    # Stash as new slot
    bcbio(x, "stash") <- "XXX"
    expect_identical(
        bcbio(x, "stash"),
        "XXX"
    )
    # Metadata stash
    bcbio(x, "metadata")[["stash"]] <- "YYY"
    expect_identical(
        bcbio(x, "metadata")[["stash"]],
        "YYY"
    )
})



# fetchPCAData =================================================================
test_that("fetchPCAData", {
    x <- fetchPCAData(pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            "pc1" = "numeric",
            "pc2" = "numeric",
            "ident" = "factor",
            "nGene" = "integer",
            "nUMI" = "integer",
            "origIdent" = "factor",
            "res0x8" = "factor",
            "res1" = "factor",
            "centerX" = "numeric",
            "centerY" = "numeric"
        )
    )
    subset <- x[1L, ] %>%
        tibble::rownames_to_column(.) %>%
        dplyr::mutate_if(
            is.numeric,
            dplyr::funs(round(., digits = 3L))
        ) %>%
        tibble::column_to_rownames(.)
    target <- data.frame(
        "pc1" = 0.443,
        "pc2" = -1.575,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = 47L,
        "nUMI" = 70L,
        "origIdent" = factor("SeuratProject"),
        "res0x8" = factor("0", levels = "0"),
        "res1" = factor("0", levels = c("0", "1", "2", "3")),
        "centerX" = -1.176,
        "centerY" = -1.683,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    # Use `glimpse()` to compare the rounded numbers
    expect_equal(subset, target)
})



# fetchTSNEData ================================================================
test_that("fetchTSNEData", {
    x <- fetchTSNEData(pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            "tSNE1" = "numeric",
            "tSNE2" = "numeric",
            "ident" = "factor",
            "nGene" = "integer",
            "nUMI" = "integer",
            "origIdent" = "factor",
            "res0x8" = "factor",
            "res1" = "factor",
            "centerX" = "numeric",
            "centerY" = "numeric"
        )
    )

    subset <- head(x, 1L) %>%
        tibble::rownames_to_column(.) %>%
        dplyr::mutate_if(
            is.numeric,
            dplyr::funs(round(., digits = 3L))
        ) %>%
        tibble::column_to_rownames(.)
    target <- data.frame(
        "tSNE1" = 14.422,
        "tSNE2" = 8.336,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = 47L,
        "nUMI" = 70L,
        "origIdent" = factor("SeuratProject"),
        "res0x8" = factor("0"),
        "res1" = factor("0", levels = c("0", "1", "2", "3")),
        "centerX" = -0.332,
        "centerY" = 18.76,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    expect_equal(subset, target)
})



# fetchTSNEExpressionData ======================================================
test_that("fetchTSNEExpressionData", {
    x <- fetchTSNEExpressionData(
        object = pbmc_small,
        genes = head(rownames(pbmc_small))
    )
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            "tSNE1" = "numeric",
            "tSNE2" = "numeric",
            "ident" = "factor",
            "nGene" = "integer",
            "nUMI" = "integer",
            "origIdent" = "factor",
            "res0x8" = "factor",
            "res1" = "factor",
            "centerX" = "numeric",
            "centerY" = "numeric",
            "mean" = "numeric",
            "median" = "numeric",
            "sum" = "numeric"
        )
    )

    subset <- x[1L, ] %>%
        tibble::rownames_to_column(.) %>%
        dplyr::mutate_if(is.numeric, dplyr::funs(round(., digits = 3L))) %>%
        tibble::column_to_rownames(.)
    # The round function will coerce integers to numerics. This is the rationale
    # for the `as.numeric()` usage below.
    target <- data.frame(
        "tSNE1" = 14.422,
        "tSNE2" = 8.336,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = as.numeric(47L),
        "nUMI" = as.numeric(70L),
        "origIdent" = factor("SeuratProject"),
        "res0x8" = factor("0"),
        "res1" = factor("0", levels = c("0", "1", "2", "3")),
        "centerX" = -0.332,
        "centerY" = 18.76,
        "mean" = 1.656,
        "median" = 0L,
        "sum" = 9.938,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = TRUE
    )
    expect_equal(subset, target)
})



# filterCells ==================================================================
test_that("filterCells", {
    x <- filterCells(bcb_small, minGenes = 0L)
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(
        dim(x),
        c(500L, 132L)
    )
    expect_is(
        metadata(x)[["filterParams"]],
        "list"
    )
    expect_is(
        metadata(x)[["filterCells"]],
        "character"
    )
    expect_is(
        metadata(x)[["filterGenes"]],
        "character"
    )
    expect_identical(
        metadata(x)[["subset"]],
        TRUE
    )
})

test_that("filterCells : Maximum parameters", {
    # This should return an object with the same dimensions
    x <- filterCells(
        bcb_small,
        minUMIs = 0L,
        maxUMIs = Inf,
        minGenes = 0L,
        maxGenes = Inf,
        maxMitoRatio = 1L,
        minNovelty = 0L,
        minCellsPerGene = 0L
    )
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(
        dim(x),
        dim(bcb_small)
    )
})

test_that("filterCells : Cutoff failures", {
    expect_error(
        filterCells(bcb_small, minUMIs = Inf),
        "No cells passed `minUMIs` cutoff"
    )
})

test_that("filterCells : Per sample cutoffs", {
    # Get the count of sample1 (run1_AGAGGATA)
    # We're applying no filtering to that sample
    metrics <- metrics(bcb_small)
    sample <- levels(metrics[["sampleID"]])[[1L]]
    expect_identical(sample, "multiplexed_AAAAAAAA")
    nCells <- length(which(metrics[["sampleID"]] == sample))
    x <- filterCells(
        object = bcb_small,
        minUMIs = c("multiplexed_AAAAAAAA" = 0L),
        maxUMIs = c("multiplexed_AAAAAAAA" = Inf),
        minGenes = c("multiplexed_AAAAAAAA" = 0L),
        maxGenes = c("multiplexed_AAAAAAAA" = Inf),
        maxMitoRatio = c("multiplexed_AAAAAAAA" = 0L),
        minNovelty = c("multiplexed_AAAAAAAA" = 0L)
    )
    expect_identical(ncol(x), nCells)
})



# metrics ======================================================================
test_that("metrics : bcbioSingleCell", {
    x <- metrics(bcb_small)
    expect_identical(
        lapply(x, class),
        list(
            "nCount" = "integer",
            "nUMI" = "integer",
            "nGene" = "integer",
            "nCoding" = "integer",
            "nMito" = "integer",
            "log10GenesPerUMI" = "numeric",
            "mitoRatio" = "numeric",
            "sampleID" = "factor",
            "sampleName" = "factor",
            "description" = "factor",
            "fileName" = "factor",
            "index" = "factor",
            "sequence" = "factor",
            "sampleNameAggregate" = "factor",
            "revcomp" = "factor",
            "interestingGroups" = "factor"
        )
    )
})

test_that("metrics : seurat", {
    x <- metrics(pbmc_small)
    expect_identical(
        lapply(x, class),
        list(
            "nGene" = "integer",
            "nUMI" = "integer",
            "origIdent" = "factor",
            "res0x8" = "factor",
            "res1" = "factor",
            "sampleID" = "factor",
            "sampleName" = "factor",
            "description" = "factor",
            "interestingGroups" = "factor",
            "ident" = "factor"
        )
    )
})



# selectSamples ================================================================
test_that("selectSamples : bcbioSingleCell", {
    x <- selectSamples(bcb_small, sampleName = "rep_1")
    expect_s4_class(x, "bcbioSingleCell")
    expect_true(metadata(x)[["selectSamples"]])
    expect_identical(
        dim(x),
        c(500L, 500L)
    )
    expect_identical(
        sampleData(x)[["sampleID"]],
        factor("multiplexed_AAAAAAAA")
    )
    # Ensure that bcbio cellular barcodes get updated correctly
    cb <- metadata(x)[["cellularBarcodes"]]
    expect_is(cb, "list")
    ids1 <- sampleData(x) %>%
        .[["sampleID"]] %>%
        as.character()
    ids2 <- names(cb)
    expect_identical(ids1, ids2)
    # Check that tibbles are stored inside list
    expect_identical(
        lapply(cb, class) %>%
            unlist() %>%
            unique(),
        c("tbl_df", "tbl", "data.frame")
    )
})

test_that("selectSamples : Match failure", {
    expect_error(
        selectSamples(bcb_small, sampleID = "XXX"),
        "sampleID metadata column doesn't contain XXX"
    )
})



# subsetPerSample ==============================================================
test_that("subsetPerSample : bcbioSingleCell", {
    subsetPerSample(bcb_small, dir = "subsetPerSample")
    expect_identical(
        list.files("subsetPerSample"),
        "multiplexed_AAAAAAAA.rda"
    )
    load("subsetPerSample/multiplexed_AAAAAAAA.rda")
    expect_identical(
        dim(multiplexed_AAAAAAAA),
        c(500L, 500L)
    )
    unlink("subsetPerSample", recursive = TRUE)
})



# topBarcodes ==================================================================
test_that("topBarcodes : bcbioSingleCell", {
    x <- topBarcodes(bcb_small)
    expect_is(x, "data.frame")
    expect_identical(
        rownames(x)[[1L]],
        "multiplexed_AAAAAAAA_CTAGCACG_AATCGGGT"
    )
})

test_that("topBarcodes : seurat", {
    x <- topBarcodes(pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        rownames(x)[[1L]],
        "GACATTCTCCACCT"
    )
})
