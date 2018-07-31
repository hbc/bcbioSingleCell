# zinbwave will calculate normalizedValues and weights matrices
# @seealso [zinbwave::zinbwave()].
.zinbwave <- function(
    Y,
    BPPARAM = BiocParallel::SerialParam(),
    epsilon = 1e12
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
            epsilon = epsilon
        )
    }))
    stopifnot(is(zinb, "SingleCellExperiment"))
    assert_are_identical(
        names(assays(zinb)),
        c("counts", "normalizedValues", "weights")
    )
    zinb
}
