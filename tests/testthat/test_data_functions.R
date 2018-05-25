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
    x <- fetchPCAData(Seurat::pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "numeric",
            res.0.8 = "character",  # factor
            res.1 = "character",  # factor
            ident = "factor",
            orig.ident = "factor",
            sampleName = "factor",
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
        sampleID = factor("SeuratProject"),
        nGene = 47L,
        nUMI = 70L,
        res.0.8 = "0",
        res.1 = "0",
        ident = factor("0", levels = c("0", "1", "2", "3")),
        orig.ident = factor("SeuratProject"),
        sampleName = factor("SeuratProject"),
        interestingGroups = factor("SeuratProject"),
        PC1 = 0.443,
        PC2 = -1.575,
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
            res.0.8 = "character",
            res.1 = "character",
            ident = "factor",
            orig.ident = "factor",
            sampleName = "factor",
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
        sampleID = factor("SeuratProject"),
        nGene = 47L,
        nUMI = 70L,
        res.0.8 = "0",
        res.1 = "0",
        ident = factor("0", levels = c("0", "1", "2", "3")),
        orig.ident = factor("SeuratProject"),
        sampleName = factor("SeuratProject"),
        interestingGroups = factor("SeuratProject"),
        tSNE_1 = 14.422,
        tSNE_2 = 8.336,
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
            res.0.8 = "character",
            res.1 = "character",
            ident = "factor",
            orig.ident = "factor",
            sampleName = "factor",
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
        sampleID = factor("SeuratProject"),
        nGene = 47L,
        nUMI = 70L,
        res.0.8 = "0",
        res.1 = "0",
        ident = factor("0", levels = c("0", "1", "2", "3")),
        orig.ident = factor("SeuratProject"),
        sampleName = factor("SeuratProject"),
        interestingGroups = factor("SeuratProject"),
        tSNE_1 = 14.422,
        tSNE_2 = 8.336,
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
        interestingGroups(Seurat::pbmc_small),
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
    expect_error(interestingGroups(x) <- "sampleName")
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
    expect_identical(lapply(x, class), all)

    # Clean mode (factor columns only)
    x <- sampleData(seurat_small, clean = TRUE)
    expect_identical(x, clean)

    # Minimal data for other seurat objects
    x <- sampleData(Seurat::pbmc_small)
    y <- DataFrame(
        "sampleName" = factor("SeuratProject"),
        row.names = "SeuratProject"
    )
    expect_identical(x, y)
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
