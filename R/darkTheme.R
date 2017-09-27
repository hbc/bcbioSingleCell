#' Dark Theme for `ggplot` Objects
#'
#' Add this last to a ggplot2 customization chain operation.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @seealso Modified version of [Seurat::DarkTheme()].
#'
#' @return No return.
#' @export
darkTheme <- function() {
    Seurat::DarkTheme() +
        theme(
            axis.line = element_line(
                colour = "white",
                linetype = "solid",
                size = 1
            ),
            axis.ticks = element_line(
                colour = "white",
                linetype = "solid",
                size = 1
            ),
            legend.justification = "center",
            legend.position = "bottom",
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
}
