#' Coerce Object
#'
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @return New class, specified by the "`to`" argument.
#'
#' @seealso `help("coerce", "methods")`.
#'
#' @examples
#' # SingleCellExperiment to list ====
#' x <- as(bcb_small, "list")
#' class(x)
#' names(x)
#'
#' # SingleCellExperiment to seurat ====
#' x <- as(bcb_small, "seurat")
#' class(x)
#' print(x)
NULL



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-SingleCellExperiment-list
setAs(
    from = "SingleCellExperiment",
    to = "list",
    function(from) {
        flatFiles(from)
    }
)



#' @rdname coerce
#' @name coerce-SingleCellExperiment-seurat
#' @section bcbioSingleCell to seurat:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object, using [as()] object
#' coercion.
setAs(
    from = "SingleCellExperiment",
    to = "seurat",
    function(from) {
        # Require that technical replicates are aggregated
        if ("sampleNameAggregate" %in% colnames(sampleData(from))) {
            inform(paste(
                "`sampleNameAggregate` metadata column detected.",
                "Use `aggregateReplicates()` to combine technical replicates.",
                sep = "\n"
            ))
        }

        # Require filtering of cells and genes
        from <- .applyFilterCutoffs(from)

        # Prepare counts matrix with gene symbols as rownames
        rawData <- assay(from)
        g2s <- gene2symbol(from)
        assertIsGene2symbol(g2s)
        rownames <- make.names(g2s[["geneName"]], unique = TRUE)
        stopifnot(!any(duplicated(rownames)))
        names(rownames) <- g2s[["geneID"]]
        rownames(rawData) <- rownames

        # Prepare metadata data.frame
        metaData <- as.data.frame(slot(from, "colData"))

        # New seurat object
        seurat <- CreateSeuratObject(
            raw.data = rawData,
            project = "bcbioSingleCell",
            # Already applied filtering cutoffs for cells and genes
            min.cells = 0L,
            min.genes = 0L,
            # Default for UMI datasets
            is.expr = 0L,
            meta.data = metaData
        )

        # Check that the dimensions match exactly
        assert_are_identical(
            x = dim(from),
            y = dim(slot(seurat, "raw.data"))
        )

        # Stash metadata and rowRanges into `misc` slot
        bcbio <- list(
            "rowRanges" = rowRanges(from),
            "metadata" = metadata(from)
        )
        slot(seurat, "misc")[["bcbio"]] <- bcbio

        seurat
    }
)
