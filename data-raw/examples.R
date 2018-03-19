library(devtools)
load_all()

# bcb_small ====================================================================
bcb_small <- loadSingleCell(
    uploadDir = "inst/extdata/harvard_indrop_v3",
    sampleMetadataFile = "inst/extdata/harvard_indrop_v3.csv",
    organism = "Homo sapiens"
)

# Subset to include only the top genes and cells
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

# Apply example filtering cutoffs
bcb_small <- filterCells(
    bcb_small,
    minUMIs = 0,
    minGenes = 0,
    maxMitoRatio = Inf,
    minNovelty = 0
)

# seurat_small =================================================================
dimsUse <- seq_len(20L)
seurat_small <- as(bcb_small, "seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = dimsUse) %>%
    RunTSNE(dims.use = dimsUse, do.fast = TRUE)

# all_markers_small ============================================================
all_markers_small <- FindAllMarkers(seurat_small)
all_markers_small_sanitized <- sanitizeMarkers(
    object = seurat_small,
    markers = all_markers_small
)
known_markers_detected <- knownMarkersDetected(
    all = all_markers_sanitized,
    known = cellTypeMarkers[["homoSapiens"]]
)

use_data(
    bcb_small,
    seurat_small,
    all_markers,
    known_markers_detected,
    compress = "xz",
    overwrite = TRUE
)
