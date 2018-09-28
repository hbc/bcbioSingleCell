# FIXME Add a way to force recalculation.
# Previously we added `metrics(recalculate = TRUE)` but a different approach
# is needed.



.calculateMetrics <- function(
    counts,
    rowRanges = NULL,
    prefilter = FALSE
) {
    assert_is_any_of(counts, c("matrix", "sparseMatrix"))
    assert_has_rows(counts)
    assert_is_any_of(rowRanges, c("GRanges", "NULL"))
    assert_is_a_bool(prefilter)

    message("Calculating cellular barcode metrics...")
    message(paste(ncol(counts), "cells detected."))

    codingGenes <- character()
    mitoGenes <- character()

    rowRangesMessage <- paste(
        "`rowRanges` is required to calculate:",
        "nCoding, nMito, mitoRatio"
    )
    missingBiotype <- function() {
        message(paste(
            "Calculating metrics without biotype information.",
            rowRangesMessage,
            sep = "\n"
        ))
    }

    # Calculate nCoding and nMito, which requires annotations.
    if (length(rowRanges) > 0L) {
        assert_is_all_of(rowRanges, "GRanges")

        setdiff <- setdiff(rownames(counts), names(rowRanges))
        if (has_length(setdiff)) {
            warning(paste(
                "Genes missing in rowRanges.",
                rowRangesMessage,
                printString(setdiff),
                sep = "\n"
            ))
            # Slot the rowRanges with empty ranges for these genes.
            # The same approach is used in `makeSummarizedExperiment()`.
            rowRanges <- suppressWarnings(c(
                rowRanges,
                emptyRanges(
                    names = setdiff,
                    mcolsNames = colnames(mcols(rowRanges))
                )
            ))
        }

        # Subset ranges to match matrix.
        assert_is_subset(rownames(counts), names(rowRanges))
        rowRanges <- rowRanges[rownames(counts)]
        rowData <- as(rowRanges, "tbl_df")
        if ("broadClass" %in% colnames(rowData)) {
            # Drop rows with NA broad class
            rowData <- filter(rowData, !is.na(!!sym("broadClass")))
            # Coding genes
            codingGenes <- rowData %>%
                filter(!!sym("broadClass") == "coding") %>%
                pull("rowname")
            message(paste(length(codingGenes), "coding genes."))
            # Mitochondrial genes
            mitoGenes <- rowData %>%
                filter(!!sym("broadClass") == "mito") %>%
                pull("rowname")
            message(paste(length(mitoGenes), "mitochondrial genes."))
        } else {
            missingBiotype()
        }
    } else {
        missingBiotype()
    }

    # Following the Seurat `seurat@meta.data` naming conventions.
    data <- tibble(
        rowname = colnames(counts),
        nUMI = colSums(counts),
        nGene = colSums(counts > 0L),
        nCoding = colSums(counts[codingGenes, , drop = FALSE]),
        nMito = colSums(counts[mitoGenes, , drop = FALSE])
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
        data <- data %>%
            filter(!is.na(UQ(sym("log10GenesPerUMI")))) %>%
            filter(!!sym("nUMI") > 0L) %>%
            filter(!!sym("nGene") > 0L)
        message(paste(
            nrow(data), "/", ncol(counts),
            "cellular barcodes passed pre-filtering",
            paste0("(", percent(nrow(data) / ncol(counts)), ")")
        ))
    }

    # Return as DataFrame, containing cell IDs in the rownames.
    data %>%
        # Enforce count columns as integers (e.g. `nUMI`).
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer) %>%
        # Coerce character vectors to factors, and drop levels.
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        as("DataFrame")
}
