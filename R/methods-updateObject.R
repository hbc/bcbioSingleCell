# if (!is.factor(metadata(object[["cell2sample"]]))) {
#     cells <- colnames(object)
#     samples <- rownames(sampleMetadata(object))
#     cell2sample <- mapCellsToSamples(
#         cells = cells,
#         samples = samples)
# }



# # Version-specific fixes
# if (metadata(object)[["version"]] == "0.0.22") {
#     if (is.data.frame(cell2sample)) {
#         cells <- as.character(cell2sample[["cellID"]])
#         samples <- as.factor(cell2sample[["sampleID"]])
#         cell2sample <- samples
#         names(cell2sample) <- cells
#     } else {
#         cell2sample <- NULL
#     }
# }
