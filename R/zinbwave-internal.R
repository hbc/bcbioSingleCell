# Check assays for zinbwave weights
.hasZinbwave <- function(object) {
    stopifnot(is(object, "SingleCellExperiment"))
    all(c("normalizedValues", "weights") %in% assayNames(object))
}



# zinbwave will calculate normalizedValues and weights matrices
# @seealso [zinbwave::zinbwave()].
.runZinbwave <- function(
    Y,  # nolint
    BPPARAM = BiocParallel::SerialParam(),  # nolint
    epsilon = 1e12,
    ...
) {
    message("Running zinbwave...")
    stopifnot(is(Y, "SingleCellExperiment"))
    Y <- as(Y, "SingleCellExperiment")  # nolint
    # zinbFit doesn't support `dgCMatrix``, so coerce counts to matrix
    assays(Y) <- list(counts = as.matrix(counts(Y)))
    message(printString(system.time({
        zinb <- zinbwave(
            Y = Y,
            K = 0L,
            BPPARAM = BPPARAM,
            epsilon = epsilon,
            ...
        )
    })))
    stopifnot(is(zinb, "SingleCellExperiment"))
    assert_are_identical(
        assayNames(zinb),
        c("counts", "normalizedValues", "weights")
    )
    zinb
}



# Stash zinbwave calculations into assays slot
.slotZinbwaveIntoAssays <- function(object, ...) {
    stopifnot(is(object, "SingleCellExperiment"))
    stopifnot(.isFiltered(object))
    zinb <- .runZinbwave(object, ...)
    assays(object)[["normalizedValues"]] <-
        assays(zinb)[["normalizedValues"]]
    assays(object)[["weights"]] <-
        assays(zinb)[["weights"]]
    object
}



# Attempt to return stashed zinbwave calcs, or recalculate
.zinbwave <- function(object, ...) {
    if (.hasZinbwave(object)) {
        message("Using stashed zinbwave weights in assays")
        counts <- counts(object)
        if (!is.matrix(counts)) {
            counts <- as.matrix(counts)
            counts(object) <- counts
        }
        object
    } else {
        .runZinbwave(object, ...)
    }
}
