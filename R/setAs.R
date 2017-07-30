# Seurat ====
# [test_seurat_object.R](https://goo.gl/GMQAdC)
setAs("bcbioSCSubset", "seurat", function(from) {
    project <- deparse(substitute(from))
    counts <- counts(from, gene2symbol = TRUE)
    rownames(counts) %>%
        head(n = 10L) %>%
        toString %>%
        paste("Gene symbols:", ., "...") %>%
        message

    # Filtering criteria
    minGenes <- metadata(from) %>%
        .[["filteringCriteria"]] %>%
        .[["minGenes"]]
    maxMitoRatio <- metadata(from) %>%
        .[["filteringCriteria"]] %>%
        .[["maxMitoRatio"]]

    object <- CreateSeuratObject(
        raw.data = counts,
        project = project,
        min.genes = minGenes)

    # Add bcbio metrics as metadata
    metrics <- metrics(from)

    # Check rowname integrity
    if (!identical(object@cell.names, rownames(metrics))) {
        stop("Metrics rowname mismatch with `cell.names` slot")
    }

    if (!identical(
        as.integer(object@meta.data[["nGene"]]),
        as.integer(metrics[["nGene"]]))) {
        stop("Gene detection mismatch")
    }
    if (!identical(
        as.integer(object@meta.data[["nUMI"]]),
        as.integer(metrics[["nUMI"]]))) {
        stop("UMI detection mismatch")
    }

    # Rows must correspond to `object@cell.names`
    object <- AddMetaData(object, metrics)
    print(glimpse(object@meta.data))

    # Normalize
    object <- NormalizeData(
        object,
        normalization.method = "LogNormalize",
        scale.factor = 10000L)

    print(object)
    object
})
