# bcbioSinglecell ====
library(bcbioSinglecell)
if (file.exists("data/run.rda")) {
    data(run)
} else {
    create_new_project()
    run <- load_run(
        upload_dir = "data/indrop_rnaseq",
        organism = "mmusculus",
        metadata = "meta/indrop_rnaseq.xlsx")
    save_data(run)
}

# knitr ====
library(knitr)
opts_chunk$set(
    audodep = TRUE,
    cache = TRUE,
    cache.lazy = FALSE,
    error = FALSE,
    fig.align = "center",
    fig.height = 7,
    fig.keep = "all",
    fig.path = "figures/",
    fig.width = 7,
    message = FALSE,
    tidy = TRUE,
    warning = FALSE)

# ggplot2 ====
library(ggplot2)
theme_set(theme_light(base_size = 14))
