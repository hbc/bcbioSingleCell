# Seurat ====
# [test_seurat_object.R](https://goo.gl/GMQAdC)
setAs("bcbioSCSubset", "seurat", function(from) {
    name <- deparse(substitute(from))
    counts <- counts(from, gene2symbol = TRUE)
    metrics <- metrics(from) %>% dotted

    # Filtering criteria
    minGenes <- metadata(from) %>%
        .[["filteringCriteria"]] %>%
        .[["minGenes"]]
    maxMitoRatio <- metadata(from) %>%
        .[["filteringCriteria"]] %>%
        .[["maxMitoRatio"]]

    new("seurat", raw.data = counts) %>%
        CreateSeuratObject(
            project = name,
            min.genes = minGenes) %>%
        AddMetaData(metrics) %>%
        SubsetData(subset.name = "mito.ratio", accept.high = maxMitoRatio)
})
