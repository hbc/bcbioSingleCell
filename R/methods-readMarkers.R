#' Read Known Markers
#'
#' @rdname readMarkers
#' @name readMarkers
#'
#' @param object Gene markers file (CSV or Excel).
#' @param gene2symbol Gene-to-symbol annotation [data.frame].
#' @param show Show [kable].
#'
#' @return [tibble].
NULL



# Methods ====
#' @rdname readMarkers
#' @export
setMethod("readMarkers", "character", function(
    object,
    gene2symbol,
    show = FALSE) {
    if (!is.data.frame(gene2symbol)) {
        stop("gene2symbol must be data.frame")
    }
    if (length(dimnames(gene2symbol)[[2L]]) != 2L) {
        stop("gene2symbol must only contain two columns")
    }
    if (!identical(
        colnames(gene2symbol),
        c("ensgene", "symbol"))) {
        stop("gene2symbol colnames must be 'ensgene', 'symbol'")
    }

    markers <- readFileByExtension(object) %>%
        camel

    # Match the markers file by Ensembl gene identifier, otherwise symbol
    if ("ensgene" %in% colnames(markers)) {
        message("Matching by gene identifier")
        markers <- markers %>%
            .[, c("cellType", "ensgene")] %>%
            .[!is.na(.[["ensgene"]]), ] %>%
            left_join(gene2symbol, by = "ensgene")
        # Check for bad identifiers
        if (any(is.na(markers[["symbol"]]))) {
            missing <- markers %>%
                .[is.na(.[["symbol"]]), ] %>%
                pull("ensgene") %>%
                sort %>%
                unique
            stop(paste("Bad genes:", toString(missing)))
        }
    } else if ("symbol" %in% colnames(markers)) {
        message("Matching by gene symbol")
        markers <- markers %>?%
            .[, c("cellType", "symbol")] %>%
            .[!is.na(.[["symbol"]]), ] %>%
            left_join(gene2symbol, by = "symbol")
        # Check for bad identifiers
        if (any(is.na(markers[["ensgene"]]))) {
            missing <- markers %>%
                .[is.na(.[["ensgene"]]), ] %>%
                pull("symbol") %>%
                sort %>%
                unique
            stop(paste("Bad symbols:", toString(missing)))
        }
    } else {
        stop("Marker file must contain 'ensgene' or 'symbol'")
    }

    markers <- markers %>%
        .[!is.na(.[["ensgene"]]), ] %>%
        .[, c("cellType", "symbol", "ensgene")] %>%
        arrange(!!!syms(c("cellType", "symbol"))) %>%
        distinct

    if (isTRUE(show)) {
        kable(markers, caption = "Known markers") %>% show
    }

    markers
})
