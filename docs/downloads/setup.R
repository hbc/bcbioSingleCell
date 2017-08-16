library(knitr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(bcbioSinglecell, quietly = TRUE)

opts_chunk$set(
    audodep = TRUE,
    cache = TRUE,
    cache.lazy = FALSE,
    error = FALSE,
    fig.height = 8,
    fig.retina = 2,
    fig.width = 8,
    message = FALSE,
    tidy = TRUE,
    warning = FALSE)

theme_set(theme_light(base_size = 14))
theme_update(legend.position = "bottom")
