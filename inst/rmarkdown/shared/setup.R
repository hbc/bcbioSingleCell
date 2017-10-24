library(bcbioSingleCell)
library(knitr)
library(viridis)
library(stringr)
library(tidyverse)

# Set seed for reproducibility
set.seed(1454944673)

opts_chunk[["set"]](
    audodep = TRUE,
    cache = FALSE,
    cache.lazy = FALSE,
    error = FALSE,
    fig.height = 10,
    fig.retina = 2,
    fig.width = 10,
    message = FALSE,
    tidy = TRUE,
    warning = FALSE)

theme_set(
    theme_light(base_size = 14))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")
