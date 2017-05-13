read_bcbio_programs <- function(run) {
    programs <- file.path(run$project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))
    return(programs)
}
