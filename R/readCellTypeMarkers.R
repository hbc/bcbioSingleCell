#' Read Cell Type Markers File
#'
#' @name readCellTypeMarkers
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param file `string`. Gene markers file path (CSV or Excel).
#'
#' @return `grouped_df`, grouped by "`cellType`" column.
#' @export
#'
#' @examples
#' # Homo sapiens
#' file <- system.file(
#'     file.path("extdata", "cell_type_markers.csv"),
#'     package = "bcbioSingleCell"
#' )
#' gene2symbol <- makeGene2symbolFromEnsembl("Homo sapiens")
#' readCellTypeMarkers(file, gene2symbol = gene2symbol)
readCellTypeMarkers <- function(file, gene2symbol) {
    assertIsGene2symbol(gene2symbol)
    data <- readFileByExtension(file) %>%
        camel() %>%
        .[, c("cellType", "geneID")] %>%
        .[complete.cases(.), , drop = FALSE]

    # Warn user about markers that aren't present in the gene2symbol.
    # This is useful for informing about putative markers that aren't expressed.
    setdiff <- setdiff(data[["geneID"]], gene2symbol[["geneID"]])
    if (length(setdiff)) {
        warning(paste(
            "Markers missing from gene2symbol (not expressed):",
            printString(setdiff),
            sep = "\n"
        ))
    }

    intersect <- intersect(data[["geneID"]], gene2symbol[["geneID"]])
    assert_is_non_empty(intersect)

    data %>%
        .[.[["geneID"]] %in% intersect, , drop = FALSE] %>%
        left_join(gene2symbol, by = "geneID") %>%
        unique() %>%
        group_by(!!sym("cellType")) %>%
        arrange(!!sym("geneName"), .by_group = TRUE)
}
