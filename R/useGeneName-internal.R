# Figure out which gene ID type is in use (geneID or geneName)
# Return TRUE if we're working with symbols (not recommended but used by Seurat)
.useGeneName <- function(object) {
    geneName <- as.character(gene2symbol(object)[["geneName"]])
    assert_is_non_empty(geneName)
    any(geneName %in% rownames(object))
}
