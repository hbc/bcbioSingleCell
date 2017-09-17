# bcbioSingleCell 0.0.18

- Renamed `plotClusters()` to `plotMarkers()`. Added soft deprecation.
- Added [viridis][] color support in tSNE plots and heatmaps.
- Converted `loadSingleCellRun()` and `loadCellRanger()` from S4 generics back to standard functions.
- Added t-SNE utility functions: `fetchTSNEData()`, `fetchTSNEExpressionData()`, and `plotTSNEExpressionData()`. This enable plotting of geometric mean values of desired marker genes.
- Updated NEWS to Markdown, with hyperlinks.
- Offloaded generics that would otherwise conflict with bcbioRNASeq to the basejump package.
- Improved roxygen documentation. Moved as much documentation as possible to the methods files.
- Updated `cellCycleMarkers` and `cellTypeMarkers` data. Now supports Drosophila.
- Sample IDs are now sanitized using `make.names()` instead of `camel()`. This avoids undesirable coercion of some IDs (e.g. `group1_1` into `group11`).
- Added recommended package syntax guidelines.
- lintr checks now allow implicit integers (e.g. `1` instead of `1L`).
- Added Seurat as dependency in `DESCRIPTION` file. The package now attaches Seurat automatically.
- Package no longer imports mononcle or suggests scater, scde, or scone. We're planning on adding these back in a future update, but build checks on Travis CI otherwise take too long.
- Added new `quantileHeatmap()` function.
- Improved Markdown header support across functions, where applicable.
- Improved `bcbioSCFiltered` to `seurat` coercion to slot relevant bcbio metadata.



# bcbioSingleCell 0.0.17

- Renamed package from `bcbioSinglecell` to `bcbioSingleCell`.
- Added [viridis][] color palette support to quality control plots.
- Added cell-cycle marker genes for human and mouse.
- Slotted `organism` in `bcbioSCDataSet` metadata, in addition to `genomeBuild`.
- Updated `bcbioSCFiltered` to `seurat` coercion method to also run `FindVariableGenes()` and `ScaleData()` by default.
- Added cell-cycle regression into [Seurat][] clustering [RMarkdown][] template.
- Improved [pkgdown][] settings and website appearance.



# bcbioSingleCell 0.0.16

- Support for CRAN release of [Seurat][].
- Improved documentation of package NAMESPACE in `bcbioSinglecell-package.R` file.
- Offloaded `download()` functionality to [basejump][] package, to avoid collisions with [bcbioRNASeq][] package. Function has been renamed to `externalFile()`.
- Removed function deprecations to simplify the NAMESPACE.
- Improved [Cell Ranger][] sample matching. With these changes the internal `.detectPipeline()` function is no longer needed.
- Improved sample directory matching in internal `.sampleDirs()` function.
- Updated paths to [Cell Ranger][] MatrixMarket files in internal `.readSparseCounts()` function.
- Updated use of `packageSE()` to `prepareSE()`, matching the corresponding [basejump][] function change.
- Renamed internal use of `filter()` to `tidy_filter()`, to avoid future NAMESPACE collisions with [ensembldb][] package.
- Renamed `bcbioSCSubset` class to `bcbioSCFiltered` class.
- Fixed memory issue in `plotZeroesVsDepth()` for datasets with high cell counts.
- Improved sample matching for `selectSamples()`. Will attempt to migrate this to bracket-based subsetting in a future update.
- Restricted sample selection with `selectSamples()` to only work on `bcbioSCFiltered` class for the time being. We can add bracket-based subsetting or S4 method support in `selectSamples()` to properly work on `bcbioSCDataSet` class in a future update.
- Updated [RMarkdown][] template settings and simplified the setup chunks.
- Initial support for `bcbioSCFiltered` class coercion to [monocle][] `CellDataSet` class.
- Suggest [scater][] and [scone][] packages for additional quality control and visualization.
- Require [monocle][] >= 2.5.0.
- Added initial support for [viridis][] color palette in quality control plots.
- Suggest [rmarkdown][] and [scde][] packages.
- Initial commit of `subsetPerSample()` function.
- Miscellaneous [RMarkdown][] template improvements.



# bcbioSingleCell 0.0.15

- Renamed functions in `lowerCamelCase` from `snake_case`.
- Draft support for [monocle][], [scater][], and [scone][].
- Package now depends on [SummarizedExperiment][].
- Renamed `loadRun()` to `loadSingleCellRun()` for improved compatibility with [bcbioRNASeq][] package. This helps avoid NAMESPACE collisions between packages.
- Improved support for custom GTF files.
- Split out [Cell Ranger][] import into a separate utility function named `loadCellRanger()`.
- Added plotZerosVsDepth() by @roryk.
- Renamed `filteringCriteria` to `filterParams` in `@metadata` slot.
- Miscellaneous documentation fixes.



# bcbioSingleCell 0.0.14

- Migrated plotting functions to S4 methods.
- Improved sample metadata handling.
- Updated sample name and cellular barcode sanization.
- Initial commit of [Seurat][] utility functions.
- Initial commit of YAML parameters for [RMarkdown][] templates.
- Draft support for [Seurat][] v2 pre-release.
- Initial support for [SureCell][] UMIs.



# bcbioSingleCell 0.0.13

- Integrated [bcbio][] and [Cell Ranger][] workflows into `load_run()`.



# bcbioSingleCell 0.0.12

- Initial support for S4 object creation using [SummarizedExperiment][].



# bcbioSingleCell 0.0.11

- Added [testthat][], [lintr][], and [covr][] support for code coverage.



# bcbioSingleCell 0.0.10

- Updated [RMarkdown][] templates.



# bcbioSingleCell 0.0.9

- Compatibility update for basejump S4 NAMESPACE changes.
- Improved plots for quality control and filtering.



# bcbioSingleCell 0.0.8

- Renamed mitochondrial plot functions.
- Changed presentation of mitochondrial abundance as ratio instead of percentage.
- Simplified NAMESPACE by offloading dependencies to basejump.



# bcbioSingleCell 0.0.7

- Initial support for loading of 10x Genomics [Cell Ranger][] output.
- Add detection of droplet method based on the metadata file.
- Draft support for import of [inDrop][] i5 index barcode counts.



# bcbioSingleCell 0.0.6

- Modified `load_run()` function for improved consistency with [bcbioRNASeq][] package.
- Initial incorporation of scater package into workflow.
- Initial commit of RMarkdown template skeletons.



# bcbioSingleCell 0.0.5

- Improvements to run object based on code in [bcbioRNASeq][] package.
- Better handling in plots for datasets with many samples.
- Labeled quality control cutoffs on the plots.



# bcbioSingleCell 0.0.4

- Prepared package for [tidyverse][] dependency updates.
- Updated HPC server detection.



# bcbioSingleCell 0.0.3

 - NAMESPACE consolidation and function cleanup.



# bcbioSingleCell 0.0.2

- Added quality control plotting functions.
- Converted all functions to standard evaluation.



# bcbioSingleCell 0.0.1

- Initial draft release.



[bcbioRNASeq]: http://bioinformatics.sph.harvard.edu/bcbioRNASeq
[Cell Ranger]: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
[inDrop]: https://1cell-bio.com
[RMarkdown]: http://rmarkdown.rstudio.com
[tidyverse]: http://www.tidyverse.org
[basejump]: http://steinbaugh.com/basejump
[covr]: https://github.com/jimhester/covr
[ensembldb]: http://bioconductor.org/packages/release/bioc/html/ensembldb.html
[lintr]: https://github.com/jimhester/lintr
[monocle]: http://cole-trapnell-lab.github.io/monocle-release/
[pkgdown]: https://github.com/hadley/pkgdown
[scater]: http://bioconductor.org/packages/release/bioc/html/scater.html
[scde]: http://bioconductor.org/packages/release/bioc/html/scde.html
[scone]: https://bioconductor.org/packages/release/bioc/html/scone.html
[Seurat]: http://satijalab.org/seurat
[SummarizedExperiment]: https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html
[SureCell]: https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/surecell-wta-ddseq.html
[testthat]: https://github.com/hadley/testthat
[viridis]: https://cran.r-project.org/web/packages/viridis/index.html
