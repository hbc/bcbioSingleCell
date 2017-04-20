#' Detect if R is running on Orchestra HPC cluster
#'
#' Utility function to enable CPU intensive tasks
#'
#' @author Michael Steinbaugh
#'
#' @export
detect_orchestra <- function() {
    if (Sys.info()[["login"]] == "root" &
        Sys.info()[["sysname"]] == "Linux" &
        any(
            Sys.getenv("CDC_JOINED_DOMAIN") == "med.harvard.edu",
            Sys.getenv("LSB_EXEC_CLUSTER") == "hms_orchestra",
            grepl("\\.orchestra$", Sys.getenv("HOSTNAME")),
            grepl("\\.orchestra$", Sys.getenv("LSB_HOSTS")),
            grepl("@MED\\.HARVARD\\.EDU$", Sys.getenv("USER_PRINCIPAL_NAME"))
        )) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
