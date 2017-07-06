# FIXME Set methods matching bcbioSCDataSet and SummarizedExperiment

#' Get gene symbols from Ensembl identifier rownames
#'
#' @rdname gene2symbol
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param object Object containing SummarizedExperiment.
#'
#' @return Ensembl gene symbol (`external_gene_name`) character vector.
.gene2symbol <- function(object) {
    df <- rowData(object) %>%
        as.data.frame %>%
        set_rownames(rownames(object)) %>%
        rownames_to_column %>%
        tidy_select(c("rowname", "ensgene", "symbol"))
    ok_symbol <- filter(df, !is.na(.data[["symbol"]]))
    na_symbol <- filter(df, is.na(.data[["symbol"]])) %>%
        mutate(ensgene = .data[["rowname"]],
               symbol = .data[["rowname"]])
    g2s <- bind_rows(ok_symbol, na_symbol) %>%
        mutate(symbol = make.unique(.data[["symbol"]])) %>%
        column_to_rownames %>%
        # Make sure the rownames order match the original object
        .[rownames(object), ]

    # Check for failures
    if (!identical(rownames(object), rownames(g2s))) {
        stop("Unexpected rowname mismatch")
    }
    if (any(duplicated(g2s[["symbol"]]))) {
        stop("gene symbols are not unique")
    }

    g2s[["symbol"]]
}
