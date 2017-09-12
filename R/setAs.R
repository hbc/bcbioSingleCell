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



# Seurat: seurat class ====
setAs("bcbioSCFiltered", "seurat", function(from) {
    project <- deparse(substitute(from))
    counts <- counts(from, gene2symbol = TRUE)

    seurat <- CreateSeuratObject(
        raw.data = counts,
        project = project,
        min.cells = 0L,
        # We already prefiltered by gene level
        min.genes = 0L,
        meta.data = metrics(from))

    # Integrity checks
    if (!identical(dim(counts), dim(seurat@data))) {
        stop("Dimension mismatch between input counts and seurat object")
    }

    # Complete the initalization steps
    seurat <- seurat %>%
        NormalizeData %>%
        FindVariableGenes(do.plot = FALSE) %>%
        ScaleData

    # Stash useful bcbio run metadata into `misc` slot
    seurat@misc[["bcbio"]] <- .fromMetadata(from)

    seurat
})
