#' Dark Theme for `ggplot` Objects
#'
#' Add this last to a ggplot2 customization chain operation.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom ggplot2 element_blank element_line element_rect margin theme
#'
#' @param ... Passthrough arguments to [ggplot2::theme()].
#'
#' @seealso Modified version of [Seurat::DarkTheme()].
#'
#' @seealso
#' - [ggplot2::theme()].
#' - [Seurat::DarkTheme()].
#'
#' @return [ggplot] theme.
#' @export
darkTheme <- function(...) {
    blackBackground <- element_rect(fill = "black")
    blackBackgroundNoBorder <- element_rect(fill = "black", size = 0)
    fontMargin <- 4
    whiteText <- element_text(
        colour = "white",
        margin = margin(
            t = fontMargin,
            r = fontMargin,
            b = fontMargin,
            l = fontMargin))
    whiteLine <- element_line(
        color = "white",
        linetype = "solid",
        size = 1)
    theme(
        axis.line = whiteLine,
        axis.text = whiteText,
        axis.ticks = whiteLine,
        axis.title = whiteText,
        legend.background = blackBackground,
        legend.box.background = blackBackgroundNoBorder,
        legend.justification = "center",
        legend.key = blackBackgroundNoBorder,
        legend.position = "bottom",
        legend.text = whiteText,
        legend.title = whiteText,
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.background = blackBackground,
        plot.subtitle = whiteText,
        plot.title = whiteText,
        complete = TRUE,
        validate = TRUE,
        ...)
}
