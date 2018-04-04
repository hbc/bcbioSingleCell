#' Read Cell Type Markers File
#'
#' @name readCellTypeMarkers
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param file Gene markers file (CSV or Excel).
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
        camel()

    # Require matching by Ensembl gene ID, not symbol
    assert_are_intersecting_sets(
        x = c("cellType", "geneID"),
        y = colnames(data)
    )

    data <- data[, c("cellType", "geneID")] %>%
        .[complete.cases(.), , drop = FALSE]

    assert_is_subset(data[["geneID"]], gene2symbol[["geneID"]])

    left_join(data, gene2symbol, by = "geneID") %>%
        as_tibble() %>%
        unique() %>%
        group_by(!!sym("cellType")) %>%
        arrange(!!sym("geneName"), .by_group = TRUE)
}
