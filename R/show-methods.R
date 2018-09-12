#' Show an Object
#'
#' @name show
#' @family S4 Object
#' @author Michael Steinbuagh
#'
#' @inherit methods::show
#'
#' @examples
#' show(indrops_small)
NULL



.show <- function(object) {
    validObject(object)

    # Extend the SingleCellExperiment method
    sce <- as(object, "SingleCellExperiment")

    return <- c(
        bold(paste(class(object), metadata(object)[["version"]])),
        "http://bioinformatics.sph.harvard.edu/bcbioSingleCell",
        "citation(\"bcbioSingleCell\")",
        separatorBar,
        paste(
            bold("Upload Dir:"),
            deparse(metadata(object)[["uploadDir"]])
        )
    )

    # FIXME CellRanger doesn't currently store `runDate`.
    if (!is.null(metadata(object)[["runDate"]])) {
        return <- c(
            return,
            paste(
                bold("Upload Date:"),
                metadata(object)[["runDate"]]
            )
        )
    }

    return <- c(
        return,
        paste(
            bold("R Load Date:"),
            metadata(object)[["date"]]
        ),
        paste(
            bold("Level:"),
            deparse(metadata(object)[["level"]])
        ),
        paste(
            bold("Organism:"),
            deparse(metadata(object)[["organism"]])
        ),
        paste(
            bold("Interesting Groups:"),
            deparse(metadata(object)[["interestingGroups"]])
        )
    )

    # sampleMetadataFile
    sampleMetadataFile <- metadata(object)[["sampleMetadataFile"]]
    if (length(sampleMetadataFile)) {
        return <- c(
            return,
            paste(bold("Metadata File:"), deparse(sampleMetadataFile))
        )
    }

    # Gene annotations
    m <- metadata(object)[["rowRangesMetadata"]]
    if (is.data.frame(m) && length(m)) {
        annotationHub <-
            m[m[["name"]] == "id", "value", drop = TRUE]
        ensemblRelease <-
            m[m[["name"]] == "ensembl_version", "value", drop = TRUE]
        genomeBuild <-
            m[m[["name"]] == "genome_build", "value", drop = TRUE]
        return <- c(
            return,
            paste(bold("AnnotationHub:"), deparse(annotationHub)),
            paste(bold("Ensembl Release:"), deparse(ensemblRelease)),
            paste(bold("Genome Build:"), deparse(genomeBuild))
        )
    }

    # GFF File
    gffFile <- metadata(object)[["gffFile"]]
    if (length(gffFile)) {
        return <- c(
            return,
            paste(bold("GFF File:"), deparse(gffFile))
        )
    }

    # Filtered counts logical
    return <- c(
        return,
        paste(bold("Filtered:"), .isFiltered(object))
    )

    # Include SingleCellExperiment show method.
    return <- c(
        return,
        separatorBar,
        capture.output(show(sce))
    )

    cat(return, sep = "\n")
}



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("bcbioSingleCell"),
    definition = .show
)



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("CellRanger"),
    definition = .show
)
