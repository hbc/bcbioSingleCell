read_bcbio_programs <- function(run) {
    file.path(run$project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))
}
