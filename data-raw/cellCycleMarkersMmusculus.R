# Tirosh et al, 2015
# http://satijalab.org/seurat/cell_cycle_vignette.html
# Ran a grep replacement in BBEdit to convert 'GENE' to 'Gene' symbol
# Cenpu = Mlf1ip

library(devtools)
library(tidyverse)
load_all()

cellCycleMarkersMmusculus <-
    read_csv(file.path(
        "data-raw",
        "cellCycleMarkersMmusculus.csv")) %>%
    as.data.frame %>%
    set_rownames(.$symbol) %>%
    symbol2gene(organism = "mouse") %>%
    rownames_to_column("ensgene") %>%
    as_tibble

use_data(cellCycleMarkersMmusculus, compress = "xz", overwrite = TRUE)
