#' Detect if R is running on Orchestra HPC cluster
#'
#' Utility function to enable CPU intensive tasks
#'
#' @author Michael Steinbaugh
#'
#' @export
detect_orchestra <- function() {
    if (Sys.info()[["sysname"]] == "Linux" &
        Sys.info()[["login"]] == "root" &
        grepl("^clarinet", Sys.info()[["nodename"]])) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
