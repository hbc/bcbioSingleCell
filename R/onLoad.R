# Follow the same order as `Depends:` in `DESCRIPTION` file
.onLoad <- function(libname, pkgname) {
    pkgs <-
        c("SingleCellExperiment",
          "Seurat",
          "basejump")
    lapply(seq_along(pkgs), function(a) {
        if (!pkgs[[a]] %in% (.packages())) {
            attachNamespace(pkgs[[a]])
        }
    })
    invisible()
}
