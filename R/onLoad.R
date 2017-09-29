.onLoad <- function(libname, pkgname) {
    packages <-
        c("SingleCellExperiment",
          "Seurat",
          "basejump")
    lapply(seq_along(packages), function(a) {
        if (!packages[[a]] %in% (.packages())) {
            attachNamespace(packages[[a]])
        }
    })
    invisible()
}
