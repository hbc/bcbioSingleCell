library(devtools)
load_all()
cell_type_markers <- head(cellTypeMarkers[["homoSapiens"]])
write_csv(cell_type_markers, path = "inst/extdata/cell_type_markers.csv")
