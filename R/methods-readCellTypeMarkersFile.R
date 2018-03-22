#' Read Cell Type Markers File
#'
#' @name readCellTypeMarkersFile
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param object Gene markers file (CSV or Excel).
#'
#' @return `grouped_df`, grouped by `cellType` column.
#'
#' @examples
#' cellTypeMarkersFile <- system.file(
#'     file.path("extdata", "cellTypeMarkers.csv"),
#'     package = "bcbioSingleCell"
#' )
#' gene2symbol <- gene2symbol("Mus musculus")
#' readCellTypeMarkersFile(cellTypeMarkersFile, gene2symbol = gene2symbol)
NULL



# Constructors =================================================================
.readCellTypeMarkersFile <- function(object, gene2symbol) {
    assertIsGene2symbol(gene2symbol)
    markers <- readFileByExtension(object) %>%
        camel(strict = FALSE)

    # Match the markers file by Ensembl gene identifier, otherwise use name
    assert_are_intersecting_sets(
        x = c("geneID", "geneName"),
        y = colnames(markers)
    )
    if ("geneID" %in% colnames(markers)) {
        inform("Matching by gene identifier")
        markers <- markers %>%
            .[, c("cellType", "geneID")] %>%
            .[!is.na(.[["geneID"]]), , drop = FALSE] %>%
            left_join(gene2symbol, by = "geneID")
        # Check for bad identifiers
        if (any(is.na(markers[["geneName"]]))) {
            missing <- markers %>%
                .[is.na(.[["geneName"]]), "geneID", drop = TRUE] %>%
                sort() %>%
                unique()
            abort(paste("Invalid genes:", toString(missing)))
        }
    } else if ("geneName" %in% colnames(markers)) {
        inform("Matching by gene name (symbol)")
        markers <- markers %>%
            .[, c("cellType", "geneName")] %>%
            .[!is.na(.[["geneName"]]), ] %>%
            left_join(gene2symbol, by = "geneName")
        # Check for bad identifiers
        if (any(is.na(markers[["geneID"]]))) {
            missing <- markers %>%
                .[is.na(.[["geneID"]]), "geneName", drop = TRUE] %>%
                sort() %>%
                unique()
            abort(paste("Invalid genes:", toString(missing)))
        }
    }

    markers %>%
        as_tibble() %>%
        .[!is.na(.[["geneID"]]), , drop = FALSE] %>%
        .[, c("cellType", "geneName", "geneID")] %>%
        unique() %>%
        group_by(!!sym("cellType")) %>%
        arrange(!!sym("geneName"), .by_group = TRUE)
}



# Methods ======================================================================
#' @rdname readCellTypeMarkersFile
#' @export
setMethod(
    "readCellTypeMarkersFile",
    signature("character"),
    .readCellTypeMarkersFile
)
