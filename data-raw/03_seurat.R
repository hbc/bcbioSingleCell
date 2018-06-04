# seurat_small
# 2018-06-04
library(devtools)
library(Seurat)
load_all()

# seurat_small =================================================================
dims_use <- seq_len(20L)
seurat_small <- cellranger_small %>%
    as("seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = dims_use) %>%
    RunTSNE(
        reduction.use = "pca",
        dims.use = dims_use,
        tsne.method = "Rtsne"
    ) %>%
    # Requires `umap-learn` python package
    RunUMAP(
        reduction.use = "pca",
        dims.use = dims_use,
        min_dist = 0.75
    )

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
