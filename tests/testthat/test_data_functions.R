context("Data Functions")



# aggregateReplicates ==========================================================
test_that("aggregateReplicates", {
    x <- aggregateReplicates(indrops_small)
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
            rowRanges = structure("GRanges", package = "GenomicRanges"),
            metadata = "list"
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
    x <- fetchPCAData(Seurat::pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "numeric",  # integer
            origIdent = "factor",
            res0x8 = "character",  # factor
            res1 = "character",  # factor
            sampleName = "factor",
            interestingGroups = "factor",
            ident = "factor",
            pc1 = "numeric",
            pc2 = "numeric",
            centerX = "numeric",
            centerY = "numeric"
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
        sampleID = factor("SeuratProject"),
        nGene = 47L,
        nUMI = 70L,
        origIdent = factor("SeuratProject"),
        res0x8 = "0",
        res1 = "0",
        sampleName = factor("SeuratProject"),
        interestingGroups = factor("SeuratProject"),
        ident = factor("0", levels = c("0", "1", "2", "3")),
        pc1 = 0.443,
        pc2 = -1.575,
        centerX = -1.176,
        centerY = -1.683,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    # Use `glimpse()` to compare the rounded numbers
    expect_equal(subset, target)
})



# fetchTSNEData ================================================================
test_that("fetchTSNEData", {
    x <- fetchTSNEData(Seurat::pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "numeric",
            origIdent = "factor",
            res0x8 = "character",
            res1 = "character",
            sampleName = "factor",
            interestingGroups = "factor",
            ident = "factor",
            tSNE1 = "numeric",
            tSNE2 = "numeric",
            centerX = "numeric",
            centerY = "numeric"
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
        sampleID = factor("SeuratProject"),
        nGene = 47L,
        nUMI = 70L,
        origIdent = factor("SeuratProject"),
        res0x8 = "0",
        res1 = "0",
        sampleName = factor("SeuratProject"),
        interestingGroups = factor("SeuratProject"),
        ident = factor("0", levels = c("0", "1", "2", "3")),
        tSNE1 = 14.422,
        tSNE2 = 8.336,
        centerX = -0.332,
        centerY = 18.76,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    expect_equal(subset, target)
})



# fetchTSNEExpressionData ======================================================
test_that("fetchTSNEExpressionData", {
    x <- fetchTSNEExpressionData(
        object = Seurat::pbmc_small,
        genes = head(rownames(Seurat::pbmc_small))
    )
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "numeric",
            origIdent = "factor",
            res0x8 = "character",
            res1 = "character",
            sampleName = "factor",
            interestingGroups = "factor",
            ident = "factor",
            tSNE1 = "numeric",
            tSNE2 = "numeric",
            centerX = "numeric",
            centerY = "numeric",
            mean = "numeric",
            median = "numeric",
            sum = "numeric"
        )
    )
    subset <- head(x, 1L) %>%
        tibble::rownames_to_column(.) %>%
        dplyr::mutate_if(is.numeric, dplyr::funs(round(., digits = 3L))) %>%
        tibble::column_to_rownames(.)
    # The round function will coerce integers to numerics. This is the rationale
    # for the `as.numeric()` usage below.
    target <- data.frame(
        sampleID = factor("SeuratProject"),
        nGene = 47L,
        nUMI = 70L,
        origIdent = factor("SeuratProject"),
        res0x8 = "0",
        res1 = "0",
        sampleName = factor("SeuratProject"),
        interestingGroups = factor("SeuratProject"),
        ident = factor("0", levels = c("0", "1", "2", "3")),
        tSNE1 = 14.422,
        tSNE2 = 8.336,
        centerX = -0.332,
        centerY = 18.76,
        mean = 0.333,
        median = 0L,
        sum = 2L,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    expect_equal(subset, target)
})



# filterCells ==================================================================
test_that("filterCells", {
    invisible(capture.output(
        x <- filterCells(
            object = indrops_small,
            minUMIs = 1000L,
            maxUMIs = Inf,
            minGenes = 100L,
            maxGenes = Inf,
            maxMitoRatio = 0.1,
            minNovelty = 0.7,
            minCellsPerGene = 3L
        )
    ))
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(dim(x), c(500L, 205L))
    expect_is(metadata(x)[["filterParams"]], "list")
    expect_is(metadata(x)[["filterCells"]], "character")
    expect_is(metadata(x)[["filterGenes"]], "character")
    expect_identical(metadata(x)[["subset"]], TRUE)
})

test_that("filterCells : Maximum parameters", {
    # This should return an object with the same dimensions
    invisible(capture.output(
        x <- filterCells(
            indrops_small,
            minUMIs = 0L,
            maxUMIs = Inf,
            minGenes = 0L,
            maxGenes = Inf,
            maxMitoRatio = 1L,
            minNovelty = 0L,
            minCellsPerGene = 0L
        )
    ))
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(dim(x), dim(indrops_small))
})

test_that("filterCells : Cutoff failures", {
    expect_error(
        filterCells(indrops_small, minUMIs = Inf),
        "No cells passed `minUMIs` cutoff"
    )
})

test_that("filterCells : Per sample cutoffs", {
    # Get the count of sample1 (run1_AGAGGATA)
    # We're applying no filtering to that sample
    sampleNames <- sampleNames(indrops_small)
    expect_identical(
        sampleNames,
        c(multiplexed_AAAAAAAA = "rep_1")
    )
    invisible(capture.output(
        x <- filterCells(
            object = indrops_small,
            minUMIs = c(rep_1 = 0L),
            maxUMIs = c(rep_1 = Inf),
            minGenes = c(rep_1 = 0L),
            maxGenes = c(rep_1 = Inf),
            maxMitoRatio = c(rep_1 = 0L),
            minNovelty = c(rep_1 = 0L)
        )
    ))
    expect_identical(
        metadata(x)[["filterParams"]],
        list(
            minUMIs = c(rep_1 = 0L),
            maxUMIs = c(rep_1 = Inf),
            minGenes = c(rep_1 = 0L),
            maxGenes = c(rep_1 = Inf),
            minNovelty = c(rep_1 = 0L),
            maxMitoRatio = c(rep_1 = 0L),
            minCellsPerGene = 10L
        )
    )
})



# metrics ======================================================================
test_that("metrics : bcbioSingleCell", {
    x <- metrics(indrops_small)
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nCount = "integer",
            nUMI = "integer",
            nGene = "integer",
            nCoding = "integer",
            nMito = "integer",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            sampleName = "factor",
            fileName = "factor",
            description = "factor",
            index = "factor",
            sequence = "factor",
            aggregate = "factor",
            revcomp = "factor",
            "interestingGroups" = "factor"
        )
    )
})

test_that("metrics : seurat", {
    x <- metrics(Seurat::pbmc_small)
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "numeric",
            origIdent = "factor",
            res0x8 = "character",
            res1 = "character",
            sampleName = "factor",
            interestingGroups = "factor",
            ident = "factor"
        )
    )
})



# selectSamples ================================================================
test_that("selectSamples : bcbioSingleCell", {
    x <- selectSamples(indrops_small, sampleName = "rep_1")
    expect_s4_class(x, "bcbioSingleCell")
    expect_true(metadata(x)[["selectSamples"]])
    expect_identical(dim(x), c(500L, 500L))
    expect_identical(
        rownames(sampleData(x)),
        "multiplexed_AAAAAAAA"
    )
})

test_that("selectSamples : Match failure", {
    expect_error(
        selectSamples(indrops_small, sampleName = "XXX"),
        "\"sampleName\" metadata column doesn't contain XXX"
    )
})



# subsetPerSample ==============================================================
test_that("subsetPerSample : bcbioSingleCell", {
    x <- subsetPerSample(indrops_small, assignAndSave = FALSE)
    expect_is(x, "list")
    expect_identical(names(x), "multiplexed_AAAAAAAA")
    subsetPerSample(
        object = indrops_small,
        assignAndSave = TRUE,
        dir = "subsetPerSample"
    )
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
test_that("topBarcodes : SingleCellExperiment", {
    x <- topBarcodes(cellranger_small)
    expect_is(x, "grouped_df")
    expect_identical(
        x[["cellID"]][[1L]],
        "aggregation_GTACGTATCTTTCCTC_1"
    )
})

test_that("topBarcodes : seurat", {
    x <- topBarcodes(Seurat::pbmc_small)
    expect_is(x, "grouped_df")
    expect_identical(
        x[["cellID"]][[1L]],
        "GACATTCTCCACCT"
    )
})
