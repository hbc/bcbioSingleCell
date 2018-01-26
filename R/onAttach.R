.onAttach <- function(libname, pkgname) {
    packages <- c(
        "bcbioBase",
        "SummarizedExperiment",
        "viridis",
        "Seurat"
    )
    lapply(packages, function(package) {
        if (!package %in% (.packages())) {
            attachNamespace(package)
        }
    })
    invisible()
}
