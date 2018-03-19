library(devtools)
load_all()

bcb_small <- loadSingleCell(
    uploadDir = "inst/extdata/harvard_indrop_v3",
    sampleMetadataFile = "inst/extdata/harvard_indrop_v3.csv",
    organism = "Homo sapiens"
)

counts <- counts(bcb_small)
countsPerGene <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- names(countsPerGene)
countsPerCell <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- names(countsPerCell)

bcb_small <- bcb_small[genes, cells]

bcb_small <- filterCells(
    agg_small,
    minUMIs = 0,
    minGenes = 0,
    maxMitoRatio = 0.25,
    minNovelty = 0.7
)

# Minimal simple Seurat working example
dimsUse <- seq_len(20L)
seurat_small <- as(filter_small, "seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = dimsUse) %>%
    RunTSNE(dims.use = dimsUse, do.fast = TRUE)

all_markers <- FindAllMarkers(seurat_small)
all_markers <- sanitizeMarkers(
    object = seurat_small,
    markers = all_markers
)
known_markers_detected <- knownMarkersDetected(
    all = all_markers,
    known = cellTypeMarkers[["hsapiens"]]
)

use_data(
    bcb_small,
    seurat_small,
    all_markers,
    known_markers_detected,
    compress = "xz",
    overwrite = TRUE
)
