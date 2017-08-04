# Seurat ====
# [test_seurat_object.R](https://goo.gl/GMQAdC)
setAs("bcbioSCSubset", "seurat", function(from) {
    project <- deparse(substitute(from))
    genomeBuild <- metadata(from)[["genomeBuild"]]
    counts <- counts(from, gene2symbol = TRUE, genomeBuild = genomeBuild)

    # Filtering criteria
    minGenes <- metadata(from) %>%
        .[["filteringCriteria"]] %>%
        .[["minGenes"]]
    maxMitoRatio <- metadata(from) %>%
        .[["filteringCriteria"]] %>%
        .[["maxMitoRatio"]]

    object <- Seurat::CreateSeuratObject(
        raw.data = counts,
        project = project,
        # We already prefiltered by gene level
        min.genes = 0L)

    # Add bcbio metrics as metadata
    metrics <- metrics(from)

    # Integrity checks
    if (!identical(dim(counts), dim(object@data))) {
        stop("Dimension mismatch between bcbioSCSubset and Seurat object")
    }
    if (!identical(colnames(counts), rownames(metrics))) {
        stop("Count matrix and metrics mismatch")
    }
    if (!identical(object@cell.names, rownames(metrics))) {
        stop("Metrics rowname mismatch with `cell.names` slot")
    }
    if (!identical(
        as.integer(object@meta.data[["nGene"]]),
        as.integer(metrics[["nGene"]]))) {
        stop("Gene detection mismatch", call. = FALSE)
    }
    if (!identical(
        as.integer(object@meta.data[["nUMI"]]),
        as.integer(metrics[["nUMI"]]))) {
        stop("UMI detection mismatch", call. = FALSE)
    }

    # Rows must correspond to `object@cell.names`
    object <- Seurat::AddMetaData(object, metrics)
    colnames(object@meta.data) %>%
        toString %>%
        paste("Seurat metadata:", .) %>%
        message

    # Normalize
    object <- Seurat::NormalizeData(
        object,
        normalization.method = "LogNormalize",
        scale.factor = 10000L)

    print(object)
    object
})
