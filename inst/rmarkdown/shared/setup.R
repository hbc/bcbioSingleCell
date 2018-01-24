suppressPackageStartupMessages(library(bcbioSingleCell))
library(knitr)
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
    message = FALSE,
    tidy = TRUE,
    warning = TRUE)

theme_set(
    theme_light(base_size = 14L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")
