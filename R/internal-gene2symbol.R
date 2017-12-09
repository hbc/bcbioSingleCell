.gene2symbol <- function(object, genes) {
    gene2symbol <- bcbio(object, "gene2symbol")
    if (is.null(gene2symbol)) {
        stop("Object doesn't contain Ensembl gene to symbol mappings")
    }
    gene2symbol %>%
        .[which(.[["ensgene"]] %in% genes), "symbol", drop = TRUE]
}
