#' Read Data Versions
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param projectDir Project directory.
#'
#' @return [tibble].
.dataVersions <- function(projectDir) {
    file <- file.path(projectDir, "data_versions.csv")
    if (!file.exists(file)) {
        warning(paste(basename(file), "missing"))
        return(NULL)
    }
    read_csv(file)
}



#' Read Log File
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param file Log file.
#'
#' @return [character].
.logFile <- function(file) {
    if (!file.exists(file)) {
        stop(paste(basename(file), "missing"))
    }
    read_lines(file)
}



#' Read Program Versions
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param projectDir Project directory.
#'
#' @return [tibble].
.programs <- function(projectDir) {
    file <- file.path(projectDir, "programs.txt")
    if (!file.exists(file)) {
        stop(paste(basename(file), "missing"))
    }
    read_delim(file, col_names = c("program", "version"), delim = ",")
}
