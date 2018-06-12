context("Data Functions")



# aggregateReplicates ==========================================================
test_that("aggregateReplicates", {
    x <- aggregateReplicates(indrops_small)
    expect_s4_class(x, "SingleCellExperiment")
    expect_identical(
        dim(x),
        c(500L, 500L)
    )
    map <- metadata(x)[["aggregateReplicates"]]
    expect_is(map, "factor")
    expect_identical(length(map), 500L)
    expect_identical(length(levels(map)), 500L)
})



# fetchPCAData =================================================================
test_that("fetchPCAData", {
    x <- fetchPCAData(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "integer",
            nCoding = "integer",
            nMito = "integer",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            orig.ident = "factor",
            res.0.8 = "character",  # factor
            ident = "factor",
            sampleName = "factor",
            index = "factor",
            interestingGroups = "factor",
            PC1 = "numeric",
            PC2 = "numeric",
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
        sampleID = factor("pbmc4k_1"),
        nGene = 2785L,
        nUMI = 12147L,
        nCoding = 11413L,
        nMito = 394L,
        log10GenesPerUMI = 0.843,
        mitoRatio = 0.032,
        orig.ident = factor("pbmc4k"),
        res.0.8 = "3",
        ident = factor("3", levels = c("0", "1", "2", "3", "4", "5")),
        sampleName = factor("pbmc4k"),
        index = factor("1"),
        interestingGroups = factor("pbmc4k"),
        PC1 = 2.44,
        PC2 = -2.471,
        centerX = 3.03,
        centerY = -2.471,
        row.names = "pbmc4k_1_AAACCTGCAGGCGATA",
        stringsAsFactors = FALSE
    )
    # Use `glimpse()` to compare the rounded numbers
    expect_equal(subset, target)
})



# fetchTSNEData ================================================================
test_that("fetchTSNEData", {
    x <- fetchTSNEData(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "integer",
            nCoding = "integer",
            nMito = "integer",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            orig.ident = "factor",
            res.0.8 = "character",
            ident = "factor",
            sampleName = "factor",
            index = "factor",
            interestingGroups = "factor",
            tSNE_1 = "numeric",
            tSNE_2 = "numeric",
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
        sampleID = factor("pbmc4k_1"),
        nGene = 2785L,
        nUMI = 12147L,
        nCoding = 11413L,
        nMito = 394L,
        log10GenesPerUMI = 0.843,
        mitoRatio = 0.032,
        orig.ident = factor("pbmc4k"),
        res.0.8 = "3",
        ident = factor("3", levels = c("0", "1", "2", "3", "4", "5")),
        sampleName = factor("pbmc4k"),
        index = factor("1"),
        interestingGroups = factor("pbmc4k"),
        tSNE_1 = 4.55,
        tSNE_2 = 2.186,
        centerX = 4.034,
        centerY = 1.805,
        row.names = "pbmc4k_1_AAACCTGCAGGCGATA",
        stringsAsFactors = FALSE
    )
    expect_equal(subset, target)
})



# fetchTSNEExpressionData ======================================================
test_that("fetchTSNEExpressionData", {
    x <- fetchTSNEExpressionData(
        object = seurat_small,
        genes = head(rownames(seurat_small))
    )
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "integer",
            nCoding = "integer",
            nMito = "integer",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            orig.ident = "factor",
            res.0.8 = "character",
            ident = "factor",
            sampleName = "factor",
            index = "factor",
            interestingGroups = "factor",
            tSNE_1 = "numeric",
            tSNE_2 = "numeric",
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
        sampleID = factor("pbmc4k_1"),
        nGene = 2785L,
        nUMI = 12147L,
        nCoding = 11413L,
        nMito = 394L,
        log10GenesPerUMI = 0.843,
        mitoRatio = 0.032,
        orig.ident = factor("pbmc4k"),
        res.0.8 = "3",
        ident = factor("3", levels = c("0", "1", "2", "3", "4", "5")),
        sampleName = factor("pbmc4k"),
        index = factor("1"),
        interestingGroups = factor("pbmc4k"),
        tSNE_1 = 4.55,
        tSNE_2 = 2.186,
        centerX = 4.034,
        centerY = 1.805,
        mean = 7.833,
        median = 2L,
        sum = 47L,
        row.names = "pbmc4k_1_AAACCTGCAGGCGATA",
        stringsAsFactors = FALSE
    )
    expect_equal(subset, target)
})



# gene2symbol ==================================================================
colnames <- c("geneID", "geneName")

test_that("gene2symbol : bcbioSingleCell", {
    x <- gene2symbol(indrops_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})

test_that("gene2symbol : seurat", {
    x <- gene2symbol(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})




# interestingGroups ============================================================
test_that("interestingGroups : bcbioSingleCell", {
    expect_identical(
        interestingGroups(indrops_small),
        "sampleName"
    )
})

test_that("interestingGroups<- : bcbioSingleCell", {
    error <- "The interesting groups \"XXX\" are not defined"
    expect_error(
        interestingGroups(indrops_small) <- "XXX",
        error
    )
    expect_error(
        interestingGroups(seurat_small) <- "XXX",
        error
    )
})

test_that("interestingGroups : seurat", {
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
})

test_that("interestingGroups<- : seurat", {
    interestingGroups(indrops_small) <- "sampleName"
    expect_identical(
        interestingGroups(indrops_small),
        "sampleName"
    )
    interestingGroups(seurat_small) <- "sampleName"
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    x <- Seurat::pbmc_small
    expect_error(interestingGroups(Seurat::pbmc_small) <- "sampleName")
})



# metrics ======================================================================
test_that("metrics : seurat", {
    # Check that metrics accessor data matches meta.data slot
    x <- metrics(seurat_small)
    y <- seurat_small@meta.data
    x <- x[, colnames(y)]
    expect_identical(x, y)
})



# sampleData ===================================================================
all <- list(
    "sampleName"  = "factor",
    "fileName"  = "factor",
    "description"  = "factor",
    "index" = "factor",
    "sequence" = "factor",
    "aggregate" = "factor",
    "revcomp" = "factor",
    "interestingGroups" = "factor"
)

clean <- DataFrame(
    "sampleName" = factor("rep_1"),
    row.names = factor("multiplexed_AAAAAAAA")
)

test_that("sampleData : bcbioSingleCell", {
    # Return all columns
    x <- sampleData(indrops_small, clean = FALSE)
    expect_identical(lapply(x, class), all)

    # Clean mode (factor columns only)
    x <- sampleData(indrops_small, clean = TRUE)
    expect_identical(x, clean)
})

test_that("sampleData : seurat", {
    # Return all columns
    x <- sampleData(seurat_small, clean = FALSE)
    expect_identical(
        lapply(x, class),
        list(
            sampleName = "factor",
            index = "factor",
            interestingGroups = "factor"
        )
    )

    # Clean mode (factor columns only)
    x <- sampleData(seurat_small, clean = TRUE)
    expect_identical(
        lapply(x, class),
        list(
            sampleName = "factor"
        )
    )

    # Minimal data for other seurat objects
    expect_identical(
        sampleData(Seurat::pbmc_small),
        DataFrame(
            "sampleName" = factor("SeuratProject"),
            row.names = "SeuratProject"
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
    expect_identical(dplyr::group_vars(x), "sampleID")
    expect_identical(
        lapply(x, class),
        list(
            cellID = "character",
            sampleID = "factor",
            nUMI = "integer",
            nGene = "integer",
            nCoding = "integer",
            nMito = "integer",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            sampleName = "factor",
            description = "factor",
            index = "factor",
            interestingGroups = "factor"
        )
    )
})

test_that("topBarcodes : seurat", {
    x <- topBarcodes(seurat_small)
    expect_is(x, "grouped_df")
    expect_identical(dplyr::group_vars(x), "sampleID")
    expect_identical(
        lapply(x, class),
        list(
            cellID = "character",
            sampleID = "factor",
            nGene = "integer",
            nUMI = "integer",
            nCoding = "integer",
            nMito = "integer",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            orig.ident = "factor",
            res.0.8 = "character",
            ident = "factor",
            sampleName = "factor",
            index = "factor",
            interestingGroups = "factor"
        )
    )
})
