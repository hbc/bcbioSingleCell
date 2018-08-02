#' Plot Dimensional Reduction
#'
#' - t-SNE: **t**-distributed **S**tochastic **N**eighbor **E**mbedding.
#' - PCA: **P**rincipal **C**omponent **A**nalysis.
#' - UMAP: **U**niform **M**anifold **A**pproximation and **P**rojection.
#'
#' @name plotDimensionalReduction
#' @family Clustering Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom BiocGenerics plotPCA
#'
#' @inheritParams general
#'
#' @note [plotUMAP()] requires the Python dependency `umap-learn`. We recommend
#' installing this with conda: `conda install -c conda-forge umap-learn`.
#'
#' @seealso
#' - [Seurat::DimPlot()].
#' - [UMAP GitHub repo](https://github.com/lmcinnes/umap).
#' - [Seurat Mouse Cell Atlas vignette](https://satijalab.org/seurat/mca.html).
#'
#' @return `ggplot`.
#'
#' @examples
#' # SingleCellExperiment ====
#' object <- cellranger_small
#'
#' # t-SNE
#' plotTSNE(object)
#' plotTSNE(object, pointsAsNumbers = TRUE)
#' plotTSNE(object, dark = TRUE)
#'
#' # UMAP
#' plotUMAP(object)
#'
#' # PCA
#' plotPCA(object)
NULL



# Constructors =================================================================
.plotReducedDims <- function(
    object,
    fun,  # e.g. fetchTSNEData
    axes,
    interestingGroups = "ident",
    color = getOption("bcbio.discrete.color", NULL),
    pointSize = getOption("bcbio.pointSize", 0.75),
    pointAlpha = getOption("bcbio.pointAlpha", 0.75),
    pointsAsNumbers = FALSE,
    label = getOption("bcbio.label", TRUE),
    labelSize = getOption("bcbio.labelSize", 6L),
    dark = getOption("bcbio.dark", FALSE),
    grid = getOption("bcbio.grid", FALSE),
    legend = getOption("bcbio.legend", TRUE),
    aspectRatio = getOption("bcbio.aspectRatio", 1L),
    title = NULL
) {
    .assertHasIdent(object)
    assert_is_function(fun)
    data <- fun(object)
    assert_is_data.frame(data)
    assert_is_character(axes)
    assert_is_subset(axes, colnames(data))
    assert_is_a_string(interestingGroups)
    assert_is_subset(interestingGroups, colnames(data))
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_number(pointSize)
    assert_is_a_number(pointAlpha)
    assert_is_a_bool(pointsAsNumbers)
    assert_is_a_bool(label)
    assert_is_a_number(labelSize)
    assert_is_a_bool(dark)
    assert_is_a_bool(legend)
    assertIsAStringOrNULL(title)

    if (interestingGroups == "ident") {
        # Seurat stores the ident from `FetchData()` as `object.ident`
        colorCol <- "ident"
    } else {
        colorCol <- interestingGroups
    }

    if (isTRUE(dark)) {
        theme <- theme_midnight
    } else {
        theme <- theme_paperwhite
    }

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym(axes[["x"]]),
            y = !!sym(axes[["y"]]),
            color = !!sym(colorCol)
        )
    ) +
        theme(
            aspect_ratio = aspectRatio,
            grid = grid
        )

    if (isTRUE(pointsAsNumbers)) {
        if (pointSize < 4L) pointSize <- 4L
        p <- p +
            geom_text(
                mapping = aes(
                    x = !!sym(axes[["x"]]),
                    y = !!sym(axes[["y"]]),
                    label = !!sym("ident"),
                    color = !!sym(colorCol)
                ),
                alpha = pointAlpha,
                size = pointSize,
                show.legend = legend
            )
    } else {
        p <- p +
            geom_point(
                alpha = pointAlpha,
                size = pointSize,
                show.legend = legend
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
                mapping = aes(
                    x = !!sym("centerX"),
                    y = !!sym("centerY"),
                    label = !!sym("ident")
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
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups = "ident",
        color = getOption("bcbio.discrete.color", NULL),
        pointSize = getOption("bcbio.pointSize", 0.75),
        pointAlpha = getOption("bcbio.pointAlpha", 0.75),
        pointsAsNumbers = FALSE,
        label = getOption("bcbio.label", TRUE),
        labelSize = getOption("bcbio.labelSize", 6L),
        dark = getOption("bcbio.dark", FALSE),
        grid = getOption("bcbio.grid", FALSE),
        legend = getOption("bcbio.legend", TRUE),
        aspectRatio = getOption("bcbio.aspectRatio", 1L),
        title = NULL
    ) {
        .plotReducedDims(
            object = object,
            fun = fetchTSNEData,
            axes = c(x = "tSNE_1", y = "tSNE_2"),
            interestingGroups = interestingGroups,
            color = color,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            pointsAsNumbers = pointsAsNumbers,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            legend = legend,
            aspectRatio = aspectRatio,
            title = title
        )
    }
)



#' @rdname plotDimensionalReduction
#' @export
setMethod(
    "plotTSNE",
    signature("seurat"),
    getMethod("plotTSNE", "SingleCellExperiment")
)



#' @rdname plotDimensionalReduction
#' @export
setMethod(
    "plotUMAP",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups = "ident",
        color = getOption("bcbio.discrete.color", NULL),
        pointSize = getOption("bcbio.pointSize", 0.75),
        pointAlpha = getOption("bcbio.pointAlpha", 0.75),
        pointsAsNumbers = FALSE,
        label = getOption("bcbio.label", TRUE),
        labelSize = getOption("bcbio.labelSize", 6L),
        dark = getOption("bcbio.dark", FALSE),
        grid = getOption("bcbio.grid", FALSE),
        legend = getOption("bcbio.legend", TRUE),
        aspectRatio = getOption("bcbio.aspectRatio", 1L),
        title = NULL
    ) {
        .plotReducedDims(
            object = object,
            fun = fetchUMAPData,
            axes = c(x = "UMAP1", y = "UMAP2"),
            interestingGroups = interestingGroups,
            color = color,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            pointsAsNumbers = pointsAsNumbers,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            legend = legend,
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
    getMethod("plotUMAP", "SingleCellExperiment")
)



#' @rdname plotDimensionalReduction
#' @export
setMethod(
    "plotPCA",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups = "ident",
        color = getOption("bcbio.discrete.color", NULL),
        pointSize = getOption("bcbio.pointSize", 0.75),
        pointAlpha = getOption("bcbio.pointAlpha", 0.75),
        pointsAsNumbers = FALSE,
        label = getOption("bcbio.label", TRUE),
        labelSize = getOption("bcbio.labelSize", 6L),
        dark = getOption("bcbio.dark", FALSE),
        grid = getOption("bcbio.grid", FALSE),
        legend = getOption("bcbio.legend", TRUE),
        aspectRatio = getOption("bcbio.aspectRatio", 1L),
        title = NULL
    ) {
        .plotReducedDims(
            object = object,
            fun = fetchPCAData,
            axes = c(x = "PC1", y = "PC2"),
            interestingGroups = interestingGroups,
            color = color,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            pointsAsNumbers = pointsAsNumbers,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            legend = legend,
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
    getMethod("plotPCA", "SingleCellExperiment")
)
