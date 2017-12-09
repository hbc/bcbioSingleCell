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
dimsUse <- 1:20
seurat <- as(filtered, "seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = dimsUse) %>%
    RunTSNE(dims.use = dimsUse, do.fast = TRUE)

seuratAllMarkersOriginal <- FindAllMarkers(seurat)
seuratAllMarkers <- sanitizeMarkers(
    seurat,
    markers = seuratAllMarkersOriginal)

knownMarkersDetected <- knownMarkersDetected(
    all = seuratAllMarkers,
    known = cellTypeMarkers[["hsapiens"]])

saveData(
    bcb,
    pooled,
    filtered,
    knownMarkersDetected,
    seurat,
    seuratAllMarkers,
    seuratAllMarkersOriginal,
    dir = extdataDir,
    compress = "xz",
    overwrite = TRUE)
