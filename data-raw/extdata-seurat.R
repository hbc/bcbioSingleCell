# seurat_small =================================================================
dimsUse <- seq_len(20L)
seurat_small <- as(indrop_small, "seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = dimsUse) %>%
    RunTSNE(dims.use = dimsUse, do.fast = TRUE)

# all_markers_small ============================================================
all_markers_small <- FindAllMarkers(seurat_small)
all_markers_small <- sanitizeMarkers(
    object = seurat_small,
    markers = all_markers_small
)

# known_markers_small ==========================================================
known_markers_small <- knownMarkersDetected(
    all = all_markers_small_sanitized,
    known = cellTypeMarkers[["homoSapiens"]]
)

# save =========================================================================
saveData(
    seurat_small,
    all_markers_small,
    known_markers_small,
    dir = "inst/extdata",
    compress = "xz"
)
