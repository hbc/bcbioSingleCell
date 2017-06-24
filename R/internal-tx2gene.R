#' Read transcript to gene (tx2gene) annotation file
#'
#' @rdname tx2gene
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param genome_build Genome build.
.tx2gene <- function(genome_build) {
    .annotable(genome_build, format = "tx2gene")
}
