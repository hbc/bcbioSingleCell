library(devtools)
load_all()
cell_type_markers <- head(cell_type_markers[["homoSapiens"]])
write_csv(cell_type_markers, path = "inst/extdata/cell_type_markers.csv")
