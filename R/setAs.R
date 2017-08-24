# monocle: CellDataSet class ====
setAs("bcbioSCFiltered", "CellDataSet", function(from) {
    # phenoData
    pd <- colData(from) %>%
        as.data.frame %>%
        as("AnnotatedDataFrame")
    # featureData
    fd <- rowData(from) %>%
        as.data.frame %>%
        # monocle requires `gene_short_name` column for gene symbols
        mutate(gene_short_name = .data[["symbol"]]) %>%
        set_rownames(rownames(from)) %>%
        as("AnnotatedDataFrame")
    cds <- newCellDataSet(
        cellData = assay(from),
        phenoData = pd,
        featureData = fd)

    # Check to make sure cells or genes aren't dropped. Use `as.numeric()`
    # conversion here because CellDataSet names the dimensions.
    if (!identical(as.numeric(dim(from)), as.numeric(dim(cds)))) {
        stop("Coercion to 'CellDataSet' class changed dimensions")
    }

    cds
})



# Seurat: seurat class ====
# [test_seurat_object.R](https://goo.gl/GMQAdC)
setAs("bcbioSCFiltered", "seurat", function(from) {
    project <- deparse(substitute(from))
    counts <- counts(from, gene2symbol = TRUE)

    object <- CreateSeuratObject(
        raw.data = counts,
        project = project,
        # We already prefiltered by gene level
        min.genes = 0L)

    # Add bcbio metrics as metadata
    metrics <- metrics(from)

    # Integrity checks
    if (!identical(dim(counts), dim(object@data))) {
        stop("Dimension mismatch between bcbioSCFiltered and Seurat object")
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
    object <- AddMetaData(object, metrics)

    # Normalize
    object <- NormalizeData(
        object,
        normalization.method = "LogNormalize",
        scale.factor = 10000L)

    print(object)
    object
})
