# Cell-Type Markers
# Last compiled on September 7, 2017 using Ensembl 90 (annotables)
library(devtools)
library(tidyverse)
library(readxl)

load_all()

# Store all the sheets in a list
file <- file.path("data-raw", "cellTypeMarkers.xlsx")
sheets <- excel_sheets(file)

cellTypeMarkers <- lapply(seq_along(sheets), function(a) {
    organism <- sheets[[a]]
    # Stop instead of warn on gene identifier match failure
    tryCatch(
        read_excel(file, sheet = organism) %>%
            group_by(symbol) %>%
            nest %>%
            as.data.frame %>%
            set_rownames(.$symbol) %>%
            symbol2gene(organism = organism) %>%
            rownames_to_column("ensgene") %>%
            as_tibble %>%
            unnest %>%
            select(cellClass, cellType, ensgene, symbol) %>%
            group_by(cellClass, cellType) %>%
            arrange(symbol, .by_group = TRUE),
        warning = function(w) {
            stop(w)
        })
}) %>%
    set_names(sheets)

use_data(cellTypeMarkers, compress = "xz", overwrite = TRUE)
