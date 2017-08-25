# Tirosh et al, 2015
# http://satijalab.org/seurat/cell_cycle_vignette.html

library(devtools)
library(tidyverse)
load_all()

# Human ====
# Ensembl 89
# FAM64A --> PIMREG
# HN1 --> JPT1
# MLF1IP --> CENPU
cellCycleMarkersHsapiens <-
    read_csv(file.path(
        "data-raw",
        "cellCycleMarkersHsapiens.csv")) %>%
    as.data.frame %>%
    set_rownames(.$symbol) %>%
    symbol2gene(organism = "human") %>%
    rownames_to_column("ensgene") %>%
    as_tibble

# Mouse ====
# Ran a grep replacement in BBEdit to convert human 'GENE' to 'Gene' symbol
# Ensembl 89
# Mlf1ip --> Cenpu
cellCycleMarkersMmusculus <-
    read_csv(file.path(
        "data-raw",
        "cellCycleMarkersMmusculus.csv")) %>%
    as.data.frame %>%
    set_rownames(.$symbol) %>%
    symbol2gene(organism = "mouse") %>%
    rownames_to_column("ensgene") %>%
    as_tibble

use_data(cellCycleMarkersHsapiens,
         cellCycleMarkersMmusculus,
         compress = "xz",
         overwrite = TRUE)
