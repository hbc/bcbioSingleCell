devtools::load_all()

extdataDir <- file.path("inst", "extdata")
uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
annotable <- annotable("Homo sapiens", release = 90)

bcb <- loadSingleCell(
    uploadDir = uploadDir,
    sampleMetadataFile = sampleMetadataFile,
    annotable = annotable)

pooled <- aggregateReplicates(bcb)
filtered <- filterCells(pooled)

# Minimal simple Seurat working example
seurat <- as(filtered, "seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = 1:20) %>%
    RunTSNE(dims.use = 1:20, do.fast = TRUE)

saveData(
    bcb,
    pooled,
    filtered,
    seurat,
    dir = extdataDir,
    compress = "xz",
    overwrite = TRUE)
