# snake_case versions of basejump utility functions
# https://github.com/steinbaugh/basejump/blob/master/R/kable.R
# https://github.com/steinbaugh/basejump/blob/master/R/makeNames.R
# https://github.com/steinbaugh/basejump/blob/master/R/setNames.R



kables <- function(list, captions = NULL) {
    output <- opts_knit$get("rmarkdown.pandoc.to")
    if (!is.null(output)) {
        tables <- lapply(seq_along(list), function(a) {
            kable(list[a], caption = captions[a])
        })
        return(asis_output(tables))
    } else {
        return(list)
    }
}



make_names <- function(character) {
    character %>%
        as.character %>%
        # Convert non-alphanumeric characters
        gsub("[^[:alnum:]]", ".", .) %>%
        # Combine multiple underscores
        gsub("[\\.]+", ".", .) %>%
        # Strip leading or trailing underscores
        gsub("(^\\.|\\.$)", "", .) %>%
        # Special names
        gsub("(m|nc|r)RNA", "\\1rna", .) %>%
        # Convert acronyms to mixed case
        gsub("([A-Z])([A-Z]+)", "\\1\\L\\2", ., perl = TRUE) %>%
        # Make first letter lowercase
        gsub("(^[A-Z]{1})", "\\L\\1", ., perl = TRUE) %>%
        # Convert camelCase
        gsub("([a-z0-9])([A-Z])", "\\1.\\L\\2", ., perl = TRUE) %>%
        # Ensure syntactically valid names
        make.names %>%
        # Lowercase everything
        tolower
}



make_names_snake <- function(character) {
    character %>%
        make_names %>%
        gsub("\\.", "_", .)
}



set_names_dot <- function(data) {
    set_names(data, make_names(colnames(data)))
}



set_names_snake <- function(data) {
    set_names(data, make_names_snake(colnames(data)))
}
