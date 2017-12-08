#' Cell to Sample Mappings
#'
#' @rdname cell2sample
#' @name cell2sample
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' load(system.file(
#'     file.path("inst", "extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("inst", "extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' cell2sample(bcb) %>% glimpse()
#'
#' # seurat
#' cell2sample(seurat) %>% glimpse()
NULL



# Constructors ====
#' Define Cell to Sample Mappings
#'
#' This function extracts `sampleID` from the `cellID` column using grep
#' matching.
#'
#' @importFrom dplyr bind_rows filter mutate
#' @importFrom magrittr set_colnames
#' @importFrom parallel mclapply
#' @importFrom stringr str_match
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @return Named character vector.
.cell2sample <- function(cells, samples) {
    cells <- as.character(cells)
    samples <- unique(as.character(samples))
    list <- mclapply(seq_along(samples), function(a) {
        pattern <- paste0("^(", samples[[a]], barcodePattern)
        match <- str_match(cells, pattern = pattern) %>%
            as.data.frame() %>%
            set_colnames(
                c("cellID",
                  "sampleID",
                  "cellularBarcode",
                  # Trailing number is for matching cellranger
                  "trailingNumber")
            ) %>%
            select(c("cellID", "sampleID")) %>%
            filter(!is.na(.data[["sampleID"]]))
        vec <- match[["sampleID"]]
        names(vec) <- match[["cellID"]]
        vec
    })
    cell2sample <- unlist(list)
    if (!identical(length(cells), length(cell2sample))) {
        stop("Failed to correctly match sample IDs to cellular barcodes",
             call. = FALSE)
    }
    as.factor(cell2sample)
}



# Methods ====
#' @rdname cell2sample
#' @param samples Sample identifiers.
#' @export
setMethod(
    "cell2sample",
    signature("character"),
    function(object, samples) {
        .cell2sample(
            cells = object,
            samples = samples
        )
    })



#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("bcbioSingleCell"),
    function(object) {
        cell2sample <- metadata(object)[["cell2sample"]]
        # Generate if no mappings are stashed
        if (is.null(cell2sample)) {
            cell2sample <- .cell2sample(
                cells = colnames(object),
                samples = rownames(sampleMetadata(object))
            )
            return(cell2sample)
        }
        # Version-specific fixes
        if (metadata(object)[["version"]] == "0.0.22") {
            # v0.0.22 stashed this as a data.frame instead
            if (is.data.frame(cell2sample)) {
                cells <- as.character(cell2sample[["cellID"]])
                samples <- as.factor(cell2sample[["sampleID"]])
                cell2sample <- samples
                names(cell2sample) <- cells
            } else {
                cell2sample <- NULL
            }
        }
        cell2sample <- cell2sample[colnames(object)]
        cell2sample <- droplevels(cell2sample)
        cell2sample
    })



#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    function(object) {
        cell2sample <- bcbio(object, "cell2sample")
        # Generate if no mappings are stashed
        if (is.null(cell2sample)) {
            cell2sample <- .cell2sample(
                cells = colnames(slot(object, "data")),
                samples = rownames(sampleMetadata(object))
            )
            return(cell2sample)
        }
        cell2sample
    })
