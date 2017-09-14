#' Coerce Object
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @seealso `help("coerce")`
#'
#' @examples
#' \dontrun{
#' seurat <- as(bcb, "seurat")
#' }
NULL



# Constructors ====
.fromMetadata <- function(from) {
    metadata(from) %>%
        .[c("version",
            "sampleMetadata",
            "filterParams",
            "annotable",
            "gene2symbol",
            "genomeBuild",
            "interestingGroups",
            "organism",
            "uploadDir")]
}



# Methods ====
#' @rdname coerce
#' @name coerce-bcbioSCFiltered-seurat
#'
#' @section [bcbioSCFiltered] to [seurat]:
#' Interally, this begins by calling [Seurat::CreateSeuratObject()] without any
#' additional filtering cutoffs, since we already applied them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new [seurat] class object, using [as()] object
#' coercion. Next, global-scaling normalization is applied to the raw counts
#' with [Seurat::NormalizeData()], which (1) normalizes the gene expression
#' measurements for each cell by the total expression, (2) multiplies this by a
#' scale factor (10,000 by default), and (3) log-transforms the result.
#' [Seurat::FindVariableGenes()] is then called, which calculates the average
#' expression and dispersion for each gene, places these genes into bins, and
#' then calculates a z-score for dispersion within each bin. This helps control
#' for the relationship between variability and average expression. Finally, the
#' genes are scaled and centered using the [Seurat::ScaleData()] function.
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
