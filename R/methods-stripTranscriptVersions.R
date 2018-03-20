#' Strip Transcript Versions
#'
#' @name stripTranscriptVersions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump stripTranscriptVersions
#'
#' @inheritParams general
#'
#' @return Object with the transcript versions in the rownames removed.
#' @export
setMethod(
    "stripTranscriptVersions",
    signature("dgCMatrix"),
    function(object) {
        rownames <- stripTranscriptVersions(rownames(object))
        rownames(object) <- rownames
        object
    }
)
