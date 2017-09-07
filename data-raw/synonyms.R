library(tidyverse)

genomes <- list(
    celegans = c("Invertebrates", "Caenorhabditis_elegans"),
    dmelanogaster = c("Invertebrates", "Drosophila_melanogaster"),
    hsapiens = c("Mammalia", "Homo_sapiens"),
    mmusculus = c("Mammalia", "Mus_musculus"))

synonyms <- lapply(seq_along(genomes), function(a) {
    file.path("ftp://ftp.ncbi.nih.gov",
              "gene",
              "DATA",
              "GENE_INFO",
              genomes[[a]][[1L]],
              paste0(genomes[[a]][[2L]], ".gene_info.gz")) %>%
        read_tsv %>%
        camel %>%
        filter(synonyms != "-") %>%
        rename(ncbi = geneId) %>%
        mutate(ensgene = str_extract(.data[["dbXrefs"]],
                                     "\\bENS[A-Z]+[0-9]{11}\\b")) %>%
        select(symbol, synonyms, description, ncbi, ensgene) %>%
        arrange(symbol)
}) %>%
    set_names(names(genomes))

devtools::use_data(synonyms, compress = "xz", overwrite = TRUE)
