if (!isTRUE(goalie::hasInternet())) {
    warning("No Internet connection detected.")
    return(invisible(NULL))
}
dir.create("cache", showWarnings = FALSE)
files <- "bcbioSingleCell_0.1.0.rds"
Map(
    f = function(remoteDir, file, envir) {
        destfile <- file.path("cache", file)
        if (!file.exists(destfile)) {
            utils::download.file(
                url = paste(remoteDir, file, sep = "/"),
                destfile = destfile
            )
        }
    },
    file = files,
    MoreArgs = list(
        "remoteDir" = bcbioSingleCellTestsUrl,
        "envir" = environment()
    )
)
rm(files)
