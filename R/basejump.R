#' Utility functions from basejump package
#'
#' @rdname basejump
#'
#' @author Michael Steinbaugh
#'
#' @param data Data type that supports name assignments



# export =======================================================================

## set_names ----

#' @rdname basejump
#' @description Set names in snake_case
#' @export
set_names_snake <- function(data) {
    set_names(data, make_names_snake(colnames(data)))
}

#' @rdname basejump
#' @description Set names in dot notation
#' @export
set_names_dot <- function(data) {
    set_names(data, make_names(colnames(data)))
}



# internal =====================================================================

## dna ----

comp <- function(dna) {
    dna <- toupper(dna)
    comp <- dna %>%
        # AT base pair swap
        gsub("A", "A1", .) %>%
        gsub("T", "A", .) %>%
        gsub("A1", "T", .) %>%
        # GC base pair swap
        gsub("G", "G1", .) %>%
        gsub("C", "G", .) %>%
        gsub("G1", "C", .)
    return(comp)
}

revcomp <- function(dna) {
    dna <- toupper(dna)
    comp <- comp(dna)
    revcomp <- strsplit(comp, "")[[1]] %>%
        .[order(seq_along(.), decreasing = TRUE)] %>%
        paste(., sep = "", collapse = "")
    return(revcomp)
}


## make_names ----

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


## knit ----

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
