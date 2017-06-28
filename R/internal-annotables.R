#' [Ensembl](http://www.ensembl.org/) annotations from [annotables]
#'
#' Quickly access gene annotations and transcript-to-gene (tx2gene) mappings
#' pre-compiled from [Ensembl](http://www.ensembl.org/) with the
#' [biomaRt](http://bioconductor.org/packages/release/bioc/html/biomaRt.html)
#' package. For gene annotables, the Entrez identifier is removed, to allow
#' for unique Ensembl gene identifiers.
#'
#' @rdname annotables
#' @keywords internal
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param genome_build Genome build. Consult the [annotables] documentation
#'   for a list of currently supported genomes.
#' @param format Desired table format, either `gene` or `tx2gene`.
#'
#' @return [data.frame] with unique rows per gene/transcript.
.annotable <- function(genome_build, format = "gene") {
    if (!is.character(genome_build)) {
        stop("Genome build must be a character vector")
    }
    if (!format %in% c("gene", "tx2gene")) {
        stop("Unsupported table format")
    }
    envir <- as.environment("package:annotables")

    # Remap genome build aliases
    if (genome_build == "hg19") {
        genome_build <- "grch37"
    } else if (genome_build == "hg38") {
        genome_build <- "grch38"
    } else if (genome_build == "mm10") {
        genome_build <- "grcm38"
    }

    message(paste("Using", genome_build, format, "annotable"))

    if (format == "gene") {
        annotable <- get(genome_build, envir = envir) %>%
            mutate(entrez = NULL) %>%
            distinct %>%
            mutate(broad_class = case_when(
                # Chromosome
                str_detect(.data[["chr"]],
                           regex("mito|mt", ignore_case = TRUE)) ~ "mito",
                # Biotype
                .data[["biotype"]] == "protein_coding" ~ "coding",
                .data[["biotype"]] %in%
                    c("known_ncrna",
                      "lincRNA",
                      "non_coding") ~ "noncoding",
                str_detect(.data[["biotype"]], "pseudo") ~ "pseudo",
                .data[["biotype"]] %in%
                    c("miRNA",
                      "misc_RNA",
                      "ribozyme",
                      "rRNA",
                      "scaRNA",
                      "scRNA",
                      "snoRNA",
                      "snRNA",
                      "sRNA") ~ "small",
                .data[["biotype"]] %in%
                    c("non_stop_decay",
                      "nonsense_mediated_decay") ~ "decaying",
                str_detect(.data[["biotype"]], "IG_") ~ "ig",
                str_detect(.data[["biotype"]], "TR_") ~ "tcr",
                TRUE ~ "other"))
    } else if (format == "tx2gene") {
        annotable <- paste(genome_build, "tx2gene", sep = "_") %>%
            get(envir = envir)
    }
    annotable %>%
        as.data.frame %>%
        set_rownames(.[[1L]])
}
