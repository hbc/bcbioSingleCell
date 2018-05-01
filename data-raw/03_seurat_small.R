# seurat_small
# 2018-04-30
library(devtools)
library(Seurat)
load_all()

# seurat_small =================================================================
dims_use <- seq_len(20L)
seurat_small <- indrops_small %>%
    as("seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = dims_use) %>%
    RunTSNE(dims.use = dims_use, do.fast = TRUE)

# all_markers_small ============================================================
all_markers_small <- FindAllMarkers(seurat_small)
all_markers_small <- sanitizeMarkers(
    object = seurat_small,
    markers = all_markers_small
)

# known_markers_small ==========================================================
known_markers_small <- knownMarkersDetected(
    object = all_markers_small,
    known = cellTypeMarkers[["homoSapiens"]]
)

# save =========================================================================
use_data(
    seurat_small,
    all_markers_small,
    known_markers_small,
    compress = "xz",
    overwrite = TRUE
)
