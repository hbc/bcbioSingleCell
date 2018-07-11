#' Sample Barcode Metrics
#'
#' @name metrics
#' @family Quality Control Functions
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
#' # SingleCellExperiment ====
#' x <- metrics(cellranger_small)
#' glimpse(x)
#'
#' # dgCMatrix ====
#' counts <- counts(cellranger_small)
#' class(counts)
#' x <- metrics(counts)
#' glimpse(x)
NULL



# Methods ======================================================================
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("matrix"),
    function(
        object,
        rowData = NULL,
        prefilter = FALSE
    ) {
        assert_has_rows(object)
        assert_is_a_bool(prefilter)

        message("Calculating cellular barcode metrics")
        message(paste(ncol(object), "cells detected"))

        codingGenes <- character()
        mitoGenes <- character()

        if (length(rowData)) {
            # Here we're allowing the user to pass in mismatched rowData, as
            # long as contains all rows in the object. Note that `rowData()`
            # return doesn't include them.
            if (hasRownames(rowData)) {
                rowData <- rowData[rownames(object), , drop = FALSE]
            }
            assert_are_identical(nrow(object), nrow(rowData))
            rownames(rowData) <- rownames(object)
            rowTbl <- rowData %>%
                as.data.frame() %>%
                as_tibble(rownames = "rowname")

            if ("broadClass" %in% colnames(rowTbl)) {
                # Drop rows with NA broad class
                rowTbl <- filter(rowTbl, !is.na(!!sym("broadClass")))
                # Coding genes
                codingGenes <- rowTbl %>%
                    filter(!!sym("broadClass") == "coding") %>%
                    pull("rowname")
                message(paste(length(codingGenes), "coding genes"))
                # Mitochondrial genes
                mitoGenes <- rowTbl %>%
                    filter(!!sym("broadClass") == "mito") %>%
                    pull("rowname")
                message(paste(length(mitoGenes), "mitochondrial genes"))
            }
        }

        # Following the Seurat `seurat@meta.data` conventions here
        tbl <- tibble(
            rowname = colnames(object),
            nUMI = colSums(object),
            nGene = colSums(object > 0L),
            nCoding = colSums(object[codingGenes, , drop = FALSE]),
            nMito = colSums(object[mitoGenes, , drop = FALSE])
        ) %>%
            mutate(
                log10GenesPerUMI = log10(!!sym("nGene")) / log10(!!sym("nUMI")),
                mitoRatio = !!sym("nMito") / !!sym("nUMI")
            ) %>%
            # Ensure `n`-prefixed count columns are integer.
            mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer)

        # Apply low stringency cellular barcode pre-filtering.
        # This keeps only cellular barcodes with non-zero genes.
        if (isTRUE(prefilter)) {
            tbl <- tbl %>%
                filter(!is.na(UQ(sym("log10GenesPerUMI")))) %>%
                filter(!!sym("nUMI") > 0L) %>%
                filter(!!sym("nGene") > 0L)
            message(paste(
                nrow(tbl), "/", ncol(object),
                "cellular barcodes passed pre-filtering",
                paste0("(", percent(nrow(tbl) / ncol(object)), ")")
            ))
        }

        # Return
        tbl %>%
            # Enforce count columns as integers (e.g. `nUMI`)
            mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer) %>%
            # Coerce character vectors to factors, and drop levels
            mutate_if(is.character, as.factor) %>%
            mutate_if(is.factor, droplevels) %>%
            as.data.frame() %>%
            column_to_rownames()
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
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }

        colData <- colData(object)
        # Calculate metrics on the fly, if not stashed in colData
        if (!"nUMI" %in% colnames(colData)) {
            metrics <- suppressMessages(metrics(
                object = counts(object),
                rowData = rowData(object),
                prefilter = FALSE
            ))
            # Keep only columns unique to colData
            setdiff <- setdiff(colnames(colData), colnames(metrics))
            colData <- colData[, setdiff]
            colData <- cbind(colData, metrics)
        }

        # Merge sample-level metadata, if stashed
        sampleData <- sampleData(
            object = object,
            interestingGroups = interestingGroups,
            clean = FALSE
        )
        if (!length(sampleData)) {
            colData[["sampleID"]] <- factor("unknown")
            colData[["sampleName"]] <- factor("unknown")
            colData[["interestingGroups"]] <- factor("unknown")
        } else {
            stopifnot(is(sampleData, "DataFrame"))
            sampleData[["sampleID"]] <- rownames(sampleData)
            # Keep only columns unique to colData
            setdiff <- setdiff(colnames(colData), colnames(sampleData))
            assert_is_non_empty(setdiff)
            colData <- colData[, setdiff]
            colData[["sampleID"]] <- cell2sample(object)
            colData[["cellID"]] <- rownames(colData)
            colData <- merge(
                x = colData,
                y = sampleData,
                by = "sampleID",
                all.x = TRUE
            ) %>%
                as.data.frame() %>%
                column_to_rownames("cellID") %>%
                .[rownames(colData), , drop = FALSE]
        }

        as.data.frame(colData)
    }
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(object, ...) {
        object %>%
            as("SingleCellExperiment") %>%
            metrics()
    }
)
