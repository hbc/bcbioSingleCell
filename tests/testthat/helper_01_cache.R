files <- list.files(
    path = "../../inst/extdata",
    pattern = "*.rda",
    full.names = TRUE
)
if (length(files)) {
    mapply(
        file = files,
        FUN = function(file, envir) {
            message(paste("Loading", basename(file)))
            load(file, envir = envir)
        },
        MoreArgs = list(envir = environment())
    )
}
