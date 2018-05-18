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
#' metrics(indrops_small) %>% glimpse()
#'
#' # SingleCellExperiment ====
#' metrics(cellranger_small) %>% glimpse()
#'
#' # dgCMatrix ====
#' counts <- counts(indrops_small)
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
    data <- as.data.frame(data)
    if (hasRownames(data)) {
        data <- rownames_to_column(data)
    }
    data %>%
        camel() %>%
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
    signature("matrix"),
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

        if (length(rowData)) {
            # Subset rowData to match counts matrix, if the rownames are
            # defined. Note that `rowData()` return doesn't include them.
            if (hasRownames(rowData)) {
                rowData <- rowData[rownames(object), , drop = FALSE]
            }
            assert_are_identical(nrow(object), nrow(rowData))
            rownames(rowData) <- rownames(object)

            if (!"broadClass" %in% colnames(rowData)) {
                warning("`broadClass` is not defined in rowData")
            }

            # Drop rows with NA broad class
            rowData <- rowData[!is.na(rowData[["broadClass"]]), , drop = FALSE]

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
                log10GenesPerUMI = log10(!!sym("nGene")) / log10(!!sym("nUMI")),
                mitoRatio = !!sym("nMito") / !!sym("nUMI")
            ) %>%
            # Ensure count columns are integer.
            # `colSums()` outputs as numeric.
            mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer)

        # Apply low stringency cellular barcode pre-filtering.
        # This keeps only cellular barcodes with non-zero genes.
        if (isTRUE(prefilter)) {
            data <- data %>%
                filter(!is.na(UQ(sym("log10GenesPerUMI")))) %>%
                filter(!!sym("nUMI") > 0L) %>%
                filter(!!sym("nGene") > 0L)
            message(paste(
                nrow(data), "/", ncol(object),
                "cellular barcodes passed pre-filtering",
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
    signature("dgCMatrix"),
    getMethod("metrics", "matrix")
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("dgTMatrix"),
    getMethod("metrics", "matrix")
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("SingleCellExperiment"),
    function(object, interestingGroups) {
        validObject(object)
        colData <- colData(object)

        # Early return NULL if there's no cell-level data
        if (!length(colData)) {
            warning("`colData()` is empty")
            return(NULL)
        }

        # Get the sample-level metadata
        sampleData <- sampleData(
            object = object,
            interestingGroups = interestingGroups,
            clean = FALSE
        )
        stopifnot(is(sampleData, "DataFrame"))
        sampleData[["sampleID"]] <- rownames(sampleData)

        # Keep only columns unique to colData
        setdiff <- setdiff(colnames(colData), colnames(sampleData))
        assert_is_non_empty(setdiff)
        colData <- colData[, setdiff]
        colData[["sampleID"]] <- cell2sample(object)
        colData[["cellID"]] <- rownames(colData)

        merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        ) %>%
            as.data.frame() %>%
            camel() %>%
            column_to_rownames("cellID") %>%
            .[rownames(colData), ]
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
