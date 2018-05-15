#' Plot Dimensional Reduction
#'
#' - t-SNE: **t**-distributed **S**tochastic **N**eighbor **E**mbedding.
#' - PCA: **P**rincipal **C**omponent **A**nalysis.
#' - UMAP: **U**niform **M**anifold **A**pproximation and **P**rojection.
#'
#' @note [plotUMAP()] requires the Python dependency `umap-learn`. We recommend
#' installing this with conda: `conda install -c conda-forge umap-learn`.
#' plotUMAP(seurat_small)
#'
#' @name plotDimensionalReduction
#' @family Clustering Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom BiocGenerics plotPCA
#'
#' @inheritParams general
#'
#' @seealso
#' - [Seurat::DimPlot()].
#' - [UMAP GitHub repo](https://github.com/lmcinnes/umap).
#' - [Satija Lab Mouse Cell Atlas clustering](https://satijalab.org/seurat/mca.html).
#'
#' @return `ggplot`.
#'
#' @examples
#' # seurat ====
#' # t-SNE
#' plotTSNE(seurat_small, dark = TRUE)
#' plotTSNE(seurat_small, dark = FALSE)
#'
#' # PCA
#' plotPCA(seurat_small)
#'
#' # UMAP
#' plotUMAP(seurat_small)
NULL



# Constructors =================================================================
.plotDimensionalReduction <- function(
    data,
    axes,
    interestingGroups = "ident",
    color = scale_color_hue(),
    pointsAsNumbers = FALSE,
    pointSize = 0.5,
    pointAlpha = 0.8,
    label = TRUE,
    labelSize = 6L,
    dark = TRUE,
    grid = TRUE,
    aspectRatio = 1L,
    title = NULL
) {
    assert_is_data.frame(data)
    assert_is_character(axes)
    assert_is_subset(axes, colnames(data))
    assert_is_a_string(interestingGroups)
    assert_is_subset(interestingGroups, colnames(data))
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_bool(pointsAsNumbers)
    assert_is_a_number(pointSize)
    assert_is_a_number(pointAlpha)
    assert_is_a_bool(label)
    assert_is_a_number(labelSize)
    assert_is_a_bool(dark)
    assertIsAStringOrNULL(title)

    if (interestingGroups == "ident") {
        # Seurat stores the ident from `FetchData()` as `object.ident`
        colorCol <- "ident"
    } else {
        colorCol <- interestingGroups
    }

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = axes[["x"]],
            y = axes[["y"]],
            color = colorCol
        )
    )

    # Put the dark theme call before the other ggplot aesthetics
    if (isTRUE(dark)) {
        p <- p +
            theme_midnight(
                aspect_ratio = aspectRatio,
                grid = grid
            )
    } else {
        p <- p +
            theme_paperwhite(
                aspect_ratio = aspectRatio,
                grid = grid
            )
    }

    if (isTRUE(pointsAsNumbers)) {
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = axes[["x"]],
                    y = axes[["y"]],
                    label = "ident",
                    color = colorCol
                ),
                alpha = pointAlpha,
                size = pointSize
            )
    } else {
        p <- p +
            geom_point(
                alpha = pointAlpha,
                size = pointSize
            )
    }

    if (isTRUE(label)) {
        if (isTRUE(dark)) {
            labelColor <- "white"
        } else {
            labelColor <- "black"
        }
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = "centerX",
                    y = "centerY",
                    label = "ident"
                ),
                color = labelColor,
                size = labelSize,
                fontface = "bold"
            )
    }

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Improve the axis breaks
    p <- p +
        scale_x_continuous(breaks = pretty_breaks(n = 4L)) +
        scale_y_continuous(breaks = pretty_breaks(n = 4L))

    p
}



# Methods ======================================================================
#' @rdname plotDimensionalReduction
#' @export
setMethod(
    "plotTSNE",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        color = scale_color_hue(),
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        dark = TRUE,
        grid = TRUE,
        aspectRatio = 1L,
        title = NULL
    ) {
        data <- fetchTSNEData(object)
        .plotDimensionalReduction(
            data = data,
            axes = c(x = "tSNE1", y = "tSNE2"),
            interestingGroups = interestingGroups,
            color = color,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            aspectRatio = aspectRatio,
            title = title
        )
    }
)



#' @rdname plotDimensionalReduction
#' @export
setMethod(
    "plotPCA",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        color = scale_color_hue(),
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        dark = TRUE,
        grid = TRUE,
        aspectRatio = 1L,
        title = NULL
    ) {
        data <- fetchPCAData(object)
        .plotDimensionalReduction(
            data = data,
            axes = c(x = "pc1", y = "pc2"),
            interestingGroups = interestingGroups,
            color = color,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            aspectRatio = aspectRatio,
            title = title
        )
    }
)



#' @rdname plotDimensionalReduction
#' @export
setMethod(
    "plotUMAP",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        color = scale_color_hue(),
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        dark = TRUE,
        grid = TRUE,
        aspectRatio = 1L,
        title = NULL
    ) {
        data <- fetchUMAPData(object)
        .plotDimensionalReduction(
            data = data,
            axes = c(x = "umap1", y = "umap2"),
            interestingGroups = interestingGroups,
            color = color,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            aspectRatio = aspectRatio,
            title = title
        )
    }
)
