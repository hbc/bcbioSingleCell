library(bcbioSingleCell)
library(Seurat)
library(knitr)
library(rmarkdown)
library(tidyverse)

# Set seed for reproducibility
set.seed(1454944673L)

opts_chunk[["set"]](
    audodep = TRUE,
    cache = TRUE,
    cache.lazy = FALSE,
    error = TRUE,
    fig.height = 10L,
    fig.retina = 2L,
    fig.width = 10L,
    message = TRUE,
    tidy = FALSE,
    warning = TRUE
)

theme_set(
    theme_paperwhite(
        base_size = 14L,
        legend_position = "bottom"
    )
)
