#' Read Cell Type Markers File
#'
#' @rdname readCellTypeMarkersFile
#' @name readCellTypeMarkersFile
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param object Gene markers file (CSV or Excel).
#' @param gene2symbol Gene-to-symbol annotation [data.frame].
#'
#' @return [tibble], gropued by `cell` column.
#'
#' @examples
#' cellTypeMarkersFile <- system.file(
#'     file.path("extdata", "cellTypeMarkers.csv"),
#'     package = "bcbioSingleCell")
#' gene2symbol <- annotable("Mus musculus", format = "gene2symbol")
#' readCellTypeMarkersFile(cellTypeMarkersFile, gene2symbol = gene2symbol)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase camel checkGene2symbol readFileByExtension
#' @importFrom dplyr arrange distinct left_join
#' @importFrom rlang !!! !! sym syms
#' @importFrom tibble as_tibble
.readCellTypeMarkersFile <- function(object, gene2symbol) {
    checkGene2symbol(gene2symbol)
    markers <- readFileByExtension(object) %>%
        camel(strict = FALSE)

    # Match the markers file by Ensembl gene identifier, otherwise symbol
    if ("ensgene" %in% colnames(markers)) {
        message("Matching by gene identifier")
        markers <- markers %>%
            .[, c("cell", "ensgene")] %>%
            .[!is.na(.[["ensgene"]]), , drop = FALSE] %>%
            left_join(gene2symbol, by = "ensgene")
        # Check for bad identifiers
        if (any(is.na(markers[["symbol"]]))) {
            missing <- markers %>%
                .[is.na(.[["symbol"]]), "ensgene", drop = TRUE] %>%
                sort() %>%
                unique()
            stop(paste("Bad genes:", toString(missing)))
        }
    } else if ("symbol" %in% colnames(markers)) {
        message("Matching by gene symbol")
        markers <- markers %>%
            .[, c("cell", "symbol")] %>%
            .[!is.na(.[["symbol"]]), ] %>%
            left_join(gene2symbol, by = "symbol")
        # Check for bad identifiers
        if (any(is.na(markers[["ensgene"]]))) {
            missing <- markers %>%
                .[is.na(.[["ensgene"]]), "symbol", drop = TRUE] %>%
                sort() %>%
                unique()
            stop(paste("Bad symbols:", toString(missing)))
        }
    } else {
        stop("Marker file must contain 'ensgene' or 'symbol'", call. = FALSE)
    }

    markers %>%
        as_tibble() %>%
        .[!is.na(.[["ensgene"]]), ] %>%
        .[, c("cell", "symbol", "ensgene")] %>%
        distinct() %>%
        group_by(!!sym("cell")) %>%
        arrange(!!!syms(c("cell", "symbol")))
}



# Methods ======================================================================
#' @rdname readCellTypeMarkersFile
#' @export
setMethod(
    "readCellTypeMarkersFile",
    signature("character"),
    .readCellTypeMarkersFile)
