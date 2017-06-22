.check_sparse <- function(sparse) {
    if (class(sparse)[[1L]] != "dgCMatrix") {
        stop("Sparse counts must be dgCMatrix")
    }
}
