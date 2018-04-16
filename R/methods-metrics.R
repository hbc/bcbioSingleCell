#' Sample Barcode Metrics
#'
#' @name metrics
#' @family Data Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#' @param rowData Data describing the rows of the object.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @return `data.frame` with cellular barcodes as rows.
#'
#' @examples
#' # bcbioSingleCell ====
#' metrics(bcb_small) %>% glimpse()
#'
#' # SingleCellExperiment ====
#' metrics(cellranger_small) %>% glimpse()
#'
#' # dgCMatrix ====
#' counts <- counts(bcb_small)
#' class(counts)
#' metrics(counts) %>% glimpse()
#'
#' # seurat ====
#' metrics(seurat_small) %>% glimpse()
#' metrics(Seurat::pbmc_small) %>% glimpse()
NULL



# Constructors =================================================================
# Ensure all columns are sanitized, and return as data.frame
.tidyMetrics <- function(data) {
    data <- camel(as.data.frame(data))
    if (has_rownames(data)) {
        data <- rownames_to_column(data)
    }
    data %>%
        # Enforce count columns as integers (e.g. `nUMI`)
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer) %>%
        # Coerce character vectors to factors, and drop levels
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        column_to_rownames()
}



# Methods ======================================================================
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("dgCMatrix"),
    function(
        object,
        rowData = NULL,
        prefilter = TRUE
    ) {
        assert_has_rows(object)
        assert_is_a_bool(prefilter)

        message("Calculating cellular barcode metrics")
        message(paste(ncol(object), "cells detected"))

        codingGenes <- character()
        mitoGenes <- character()

        if (!is.null(rowData)) {
            assert_are_intersecting_sets(
                x = rownames(object),
                y = rownames(rowData)
            )

            # Subset rowData to match counts matrix
            rowData <- rowData[rownames(object), , drop = FALSE]

            # Coding genes
            codingGenes <- rowData %>%
                .[.[["broadClass"]] == "coding", , drop = FALSE] %>%
                rownames() %>%
                na.omit()
            message(paste(length(codingGenes), "coding genes detected"))

            # Mitochondrial genes
            mitoGenes <- rowData %>%
                .[.[["broadClass"]] == "mito", , drop = FALSE] %>%
                rownames() %>%
                na.omit()
            message(paste(length(mitoGenes), "mitochondrial genes detected"))
        }

        data <- tibble(
            "rowname" = colnames(object),
            # Follow the Seurat `seurat@data.info` conventions
            "nUMI" = colSums(object),
            "nGene" = colSums(object > 0L),
            "nCoding" = colSums(object[codingGenes, , drop = FALSE]),
            "nMito" = colSums(object[mitoGenes, , drop = FALSE])
        ) %>%
            mutate(
                log10GenesPerUMI = log10(.data[["nGene"]]) /
                    log10(.data[["nUMI"]]),
                # Using `nUMI` here like in Seurat example
                mitoRatio = .data[["nMito"]] / .data[["nUMI"]]
            ) %>%
            # Ensure count columns are integer.
            # `colSums()` outputs as numeric.
            mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer)

        # Apply low stringency cellular barcode pre-filtering.
        # This keeps only cellular barcodes with non-zero genes.
        if (isTRUE(prefilter)) {
            data <- data %>%
                .[.[["nUMI"]] > 0L, , drop = FALSE] %>%
                .[.[["nGene"]] > 0L, , drop = FALSE] %>%
                .[!is.na(.[["log10GenesPerUMI"]]), , drop = FALSE]
            message(paste(
                nrow(data), "cellular barcodes passed pre-filtering",
                paste0("(", percent(nrow(data) / ncol(object)), ")")
            ))
        }

        .tidyMetrics(data)
    }
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("dgTMatrix"),
    getMethod("metrics", "dgCMatrix")
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("SingleCellExperiment"),
    function(object, interestingGroups) {
        sampleData <- sampleData(
            object = object,
            interestingGroups = interestingGroups
        )
        colData <- colData(object)
        # Drop any duplicate metadata columns from colData, if necessary.
        # This step was added so we can slot rich metadata in seurat objects
        # inside `object@meta.data`, otherwise this line can be removed.
        colData <- colData[setdiff(colnames(colData), colnames(sampleData))]
        assert_are_disjoint_sets(
            x = colnames(sampleData),
            y = colnames(colData)
        )

        cell2sample <- cell2sample(object)
        assert_is_factor(cell2sample)

        sampleID <- DataFrame("sampleID" = cell2sample)

        colData <- cbind(colData, sampleID)
        colData[["rowname"]] <- rownames(colData)

        data <- merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )

        # Ensure the numeric metrics columns appear first
        data <- data[, unique(c(colnames(colData), colnames(data)))]

        .tidyMetrics(data)
    }
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(object, ...) {
        fun <- getMethod("metrics", "SingleCellExperiment")
        data <- fun(object, ...)
        # Add ident column
        data[["ident"]] <- slot(object, "ident")
        data
    }
)
