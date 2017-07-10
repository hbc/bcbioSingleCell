# bcbioSinglecell ====
library(bcbioSinglecell)
if (file.exists("data/bcb.rda")) {
    data(bcb)
} else {
    bcb <- load_run(
        upload_dir = file.path("data", "singlecell"),
        sample_metadata_file = file.path("meta", "sample_barcodes.xlsx"),
        interesting_groups = "genotype")
    save_data(bcb, compress = FALSE)
}



# knitr ====
library(knitr)
opts_chunk$set(
    audodep = TRUE,
    cache = TRUE,
    cache.lazy = FALSE,
    dev = c("png", "pdf", "svg"),
    error = FALSE,
    fig.height = 8,
    fig.retina = 2,
    fig.width = 8,
    message = FALSE,
    tidy = TRUE,
    warning = FALSE)



# ggplot2 ====
library(ggplot2)
theme_set(theme_light(base_size = 14))
