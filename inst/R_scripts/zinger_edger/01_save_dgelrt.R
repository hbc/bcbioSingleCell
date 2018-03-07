# Differential Expression Per Cluster Using zingeR-edgeR
#
# Michael Steinbaugh
# 2018-03-07
#
# Compatible with bcbioSingleCell v0.0.32
# zingeR vignette v2: https://goo.gl/4rTK1w

library(parallel)
library(edgeR)
library(zingeR)
library(Seurat)
library(bcbioSingleCell)
library(tidyverse)

data_dir <- path("dir", Sys.Date())
loadData(seurat, dir = data_dir)

# Get the cell cluster identities from the final Seurat object
interestingGroups(seurat) <- "<INTGROUP>"
sampleData(seurat)
metrics(seurat) %>% glimpse()

ident <- levels(seurat@ident)
treatment <- "<TREATMENT>"
control <- "<CONTROL>"

dgelrt <- mcmapply(
    ident = ident,
    MoreArgs = list(
        object = seurat,
        treatment = treatment,
        control = control,
        prefix = "<PREFIX>",
        dir = data_dir
    ),
    FUN = function(
        object,
        ident,
        treatment,
        control,
        maxit = 1000L,
        prefix = "zinger_edger_lrt_cluster",
        envir = parent.frame(),
        dir = ".") {
        map <- metrics(object) %>%
            as_tibble() %>%
            rownames_to_column("cellID") %>%
            .[, c("ident", "interestingGroups", "cellID")] %>%
            .[.[["ident"]] == ident, ,] %>%
            group_by(!!!syms(c("ident", "interestingGroups"))) %>%
            arrange(!!sym("cellID"), .by_group = TRUE)
        numerator <- map %>%
            .[.[["interestingGroups"]] == treatment, ] %>%
            pull("cellID")
        denominator <- map %>%
            .[.[["interestingGroups"]] == control, ] %>%
            pull("cellID")
        assert_are_set_equal(map[["cellID"]], c(numerator, denominator))

        # zingeR-edgeR differential expression
        lrt <- diffExp(
            object,
            numerator = numerator,
            denominator = denominator,
            maxit = maxit)

        # Save the individual DGELRT to disk, as a backup in case parallel
        # call is suspended on cluster
        assignAndSaveData(
            name = paste(prefix, ident, sep = "_"),
            object = lrt,
            envir = envir,
            dir = dir)

        lrt
    },
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE)
saveData(dgelrt, dir = dat_dir)
