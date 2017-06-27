# knitr ====
library(knitr)
opts_chunk$set(
    audodep = TRUE,
    cache = TRUE,
    cache.lazy = FALSE,
    dev = c("png", "pdf", "svg"),
    error = FALSE,
    fig.height = 7,
    fig.retina = 2,
    fig.width = 7,
    message = FALSE,
    tidy = TRUE,
    warning = FALSE)

# ggplot2 ====
library(ggplot2)
theme_set(theme_light(base_size = 14))

# bcbioSinglecell ====
library(bcbioSinglecell)
if (file.exists("data/bcb.rda")) {
    data(bcb)
} else {
    bcb <- load_run(
        upload_dir = file.path("data", "indrop_rnaseq"),
        metadata_file = file.path("meta", "indrop_rnaseq.xlsx"),
        interesting_groups = "genotype")
    save_data(bcb, compress = FALSE)
}
