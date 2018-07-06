#' Cell-Cycle Markers
#'
#' @aliases cellCycleMarkers
#' @family Marker Data
#' @author Michael Steinbaugh
#'
#' @examples
#' names(cell_cycle_markers)
"cell_cycle_markers"



#' @rdname cell_cycle_markers
#' @usage NULL
#' @export
cellCycleMarkers <- bcbioSingleCell::cell_cycle_markers



#' Cell-Type Markers
#'
#' @aliases cellTypeMarkers
#' @family Marker Data
#' @author Michael Steinbaugh
#'
#' @examples
#' names(cell_type_markers)
"cell_type_markers"



#' @rdname cell_type_markers
#' @usage NULL
#' @export
cellTypeMarkers <- bcbioSingleCell::cell_type_markers



# Minimal Examples =============================================================
#' Sanitized Markers Example
#'
#' @family Minimal Example Data
#' @author Michael Steinbaugh
#'
#' @examples
#' glimpse(all_markers_small)
"all_markers_small"



#' Cell Ranger Example
#'
#' @family Minimal Example Data
#' @author Michael Steinbaugh
#'
#' @examples
#' show(cellranger_small)
"cellranger_small"



#' inDrops Example
#'
#' @family Minimal Example Data
#' @author Michael Steinbaugh
#'
#' @examples
#' show(indrops_small)
"indrops_small"



#' Known Markers Example
#'
#' @family Minimal Example Data
#' @author Michael Steinbaugh
#'
#' @examples
#' glimpse(known_markers_small)
"known_markers_small"



#' Seurat Example
#'
#' @family Minimal Example Data
#' @author Michael Steinbaugh
#'
#' @examples
#' show(seurat_small)
"seurat_small"
