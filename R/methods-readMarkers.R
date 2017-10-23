#' Read Known Markers
#'
#' @rdname readMarkers
#' @name readMarkers
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param object Gene markers file (CSV or Excel).
#' @param gene2symbol Gene-to-symbol annotation [data.frame].
#'
#' @return [tibble].
NULL



# Constructors ====
#' @importFrom rlang syms !!!
.readMarkers <- function(object, gene2symbol) {
    if (!is.data.frame(gene2symbol)) {
        stop("gene2symbol must be data.frame")
    }
    if (length(dimnames(gene2symbol)[[2]]) != 2) {
        stop("gene2symbol must only contain two columns")
    }
    if (!identical(
        colnames(gene2symbol),
        c("ensgene", "symbol"))) {
        stop("gene2symbol colnames must be 'ensgene', 'symbol'")
    }

    markers <- readFileByExtension(object) %>%
        camel()

    # Match the markers file by Ensembl gene identifier, otherwise symbol
    if ("ensgene" %in% colnames(markers)) {
        message("Matching by gene identifier")
        markers <- markers %>%
            .[, c("cell", "ensgene")] %>%
            .[!is.na(.[["ensgene"]]), ] %>%
            left_join(gene2symbol, by = "ensgene")
        # Check for bad identifiers
        if (any(is.na(markers[["symbol"]]))) {
            missing <- markers %>%
                .[is.na(.[["symbol"]]), ] %>%
                pull("ensgene") %>%
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
                .[is.na(.[["ensgene"]]), ] %>%
                pull("symbol") %>%
                sort() %>%
                unique()
            stop(paste("Bad symbols:", toString(missing)))
        }
    } else {
        stop("Marker file must contain 'ensgene' or 'symbol'")
    }

    markers %>%
        .[!is.na(.[["ensgene"]]), ] %>%
        .[, c("cell", "symbol", "ensgene")] %>%
        arrange(!!!syms(c("cell", "symbol"))) %>%
        distinct()
}



# Methods ====
#' @rdname readMarkers
#' @export
setMethod(
    "readMarkers",
    signature("character"),
    .readMarkers)
