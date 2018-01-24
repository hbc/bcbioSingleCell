.convertGenesToSymbols <- function(object, genes) {
    gene2symbol <- gene2symbol(object)
    if (is.null(gene2symbol)) {
        abort("NULL gene to symbol mappings detected")
    }
    symbols <- gene2symbol %>%
        .[which(.[["ensgene"]] %in% genes), "symbol", drop = TRUE]
    if (!identical(length(genes), length(symbols))) {
        abort("Failed to map all gene identifiers to symbols")
    }
    names(symbols) <- genes
    symbols
}
