# zinbwave will calculate normalizedValues and weights matrices
# @seealso [zinbwave::zinbwave()].
.zinbwave <- function(
    Y,
    BPPARAM = BiocParallel::SerialParam(),
    epsilon = 1e12,
    ...
) {
    message("Running zinbwave")
    stopifnot(is(Y, "SingleCellExperiment"))
    Y <- as(Y, "SingleCellExperiment")
    # zinbFit doesn't support `dgCMatrix``, so coerce counts to matrix
    assays(Y) <- list(counts = as.matrix(counts(Y)))
    print(system.time({
        zinb <- zinbwave(
            Y = Y,
            K = 0L,
            BPPARAM = BPPARAM,
            epsilon = epsilon,
            ...
        )
    }))
    stopifnot(is(zinb, "SingleCellExperiment"))
    assert_are_identical(
        assayNames(zinb),
        c("counts", "normalizedValues", "weights")
    )
    zinb
}



# Stash zinbwave calculations into assays slot
.zinbwaveIntoAssays <- function(object, ...) {
    stopifnot(is(object, "SingleCellExperiment"))
    stopifnot(.isFiltered(object))
    zinb <- .zinbwave(object, ...)
    assays(object)[["normalizedValues"]] <-
        assays(zinb)[["normalizedValues"]]
    assays(object)[["weights"]] <-
        assays(zinb)[["weights"]]
    object
}
