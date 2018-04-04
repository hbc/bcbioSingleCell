.tidyMetrics <- function(object) {
    assert_is_data.frame(object)
    # Ensure that rownames are set as a column before performing this chain
    object %>%
        # Enforce count columns as integers (e.g. `nUMI`)
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer) %>%
        # Coerce character vectors to factors, and drop levels
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels)
}
