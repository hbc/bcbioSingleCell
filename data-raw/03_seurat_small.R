# seurat_small
# 2018-05-14
library(devtools)
load_all()
library(Seurat)

# seurat_small =================================================================
dims_use <- seq_len(20L)
seurat_small <- indrops_small %>%
    as("seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = dims_use) %>%
    RunTSNE(dims.use = dims_use, tsne.method = "Rtsne")

# all_markers_small ============================================================
# MAST, zinbwave/DESeq2 are recommended
all_markers_small <- FindAllMarkers(seurat_small, test.use = "wilcox")
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
