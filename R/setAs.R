# Constructors ====
.fromMetadata <- function(from) {
    list(
        version = metadata(from)[["version"]],
        annotable = metadata(from)[["annotable"]],
        filterParams = metadata(from)[["filterParams"]],
        gene2symbol = metadata(from)[["gene2symbol"]],
        genomeBuild = metadata(from)[["genomeBuild"]],
        interestingGroups = metadata(from)[["interestingGroups"]],
        organism = metadata(from)[["organism"]],
        uploadDir = metadata(from)[["uploadDir"]])
}



# monocle: CellDataSet class ====
setAs("bcbioSCFiltered", "CellDataSet", function(from) {
    # CellDataSet currently extends the ExpressionSet class. We can stash
    # useful metadata in the `experimentData` slot in a future update.

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

    seurat <- CreateSeuratObject(
        raw.data = counts,
        project = project,
        # We already prefiltered by gene level
        min.genes = 0L)

    # Add bcbio metrics as metadata
    metrics <- metrics(from)

    # Integrity checks
    if (!identical(dim(counts), dim(seurat@data))) {
        stop("Dimension mismatch between bcbioSCFiltered and Seurat object")
    }
    if (!identical(colnames(counts), rownames(metrics))) {
        stop("Count matrix and metrics mismatch")
    }
    if (!identical(seurat@cell.names, rownames(metrics))) {
        stop("Metrics rowname mismatch with `cell.names` slot")
    }
    if (!identical(
        as.integer(seurat@meta.data[["nGene"]]),
        as.integer(metrics[["nGene"]]))) {
        stop("Gene detection mismatch", call. = FALSE)
    }
    if (!identical(
        as.integer(seurat@meta.data[["nUMI"]]),
        as.integer(metrics[["nUMI"]]))) {
        stop("UMI detection mismatch", call. = FALSE)
    }

    # Complete the initalization steps
    seurat <- seurat %>%
        AddMetaData(metrics) %>%
        NormalizeData %>%
        FindVariableGenes(do.plot = FALSE) %>%
        ScaleData

    # Stash useful bcbio run metadata into `misc` slot
    seurat@misc[["bcbio"]] <- .fromMetadata(from)

    seurat
})
