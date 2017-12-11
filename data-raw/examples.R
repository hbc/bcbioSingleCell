devtools::load_all()

extdataDir <- file.path("inst", "extdata")
uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
annotable <- annotable("Homo sapiens", release = 90)

bcb <- loadSingleCell(
    uploadDir = uploadDir,
    sampleMetadataFile = sampleMetadataFile,
    annotable = annotable)

# Make the example bcbioSingleCell object more minimal to save disk space
# First, let's pick only the top 250 genes that have the most robust expression.
counts <- counts(bcb)
if (!is(counts, "dgCMatrix")) {
    stop("counts should be sparse dgCMatrix")
}

nGenes <- 1000L
countsPerGene <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = nGenes)
topGenes <- names(countsPerGene) %>%
    sort()
length(topGenes)

nCells <- 500L
countsPerCell <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = nCells)
topCells <- names(countsPerCell) %>%
    sort()
length(topCells)

if (!identical(
    c(nGenes, nCells),
    c(length(topGenes), length(topCells))
)) {
    stop("Example should have same number of genes and cells")
}

bcb <- bcb[topGenes, topCells]
if (!identical(c(nGenes, nCells), dim(bcb))) {
    stop("Subset operation failed")
}
bcb
object.size(bcb) %>%
    format(units = "auto")

pooled <- aggregateReplicates(bcb)
pooled

filtered <- filterCells(
    pooled,
    maxMitoRatio = 0.25,
    minNovelty = 0.7)
filtered

# Minimal simple Seurat working example
seurat <- as(filtered, "seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(dims.use = 1:20) %>%
    RunTSNE(dims.use = 1:20, do.fast = TRUE)
seurat

saveData(
    bcb,
    pooled,
    filtered,
    seurat,
    dir = extdataDir,
    compress = "xz",
    overwrite = TRUE)
