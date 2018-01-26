.onAttach <- function(libname, pkgname) {
    imports <- c(
        "bcbioBase",
        "Seurat"
    )
    invisible(lapply(
        X = imports,
        FUN = require,
        character.only = TRUE
    ))
}
