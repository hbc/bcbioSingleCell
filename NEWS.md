# bcbioSingleCell 0.0.23

- Improved facet wrapping of aggregated samples (`sampleNameAggregate` present in sample metadata), but removing code support for wrapping by multiplexed FASTQ description.



# bcbioSingleCell 0.0.22 (2017-11-17)

- Now internally stashing a `cell2sample` data.frame, which helps speed up operations on cellular barcode metrics calculations for quality control plots.
- Improved support for optional `annotable`, `ensemblVersion`, `gtfFile`, and `sampleMetadataFile` arguments in `loadSingleCell()` function.
- Simplified some of the messages shown during sample loading, in an attempt to make them clearer and more informative. Also improved messages shown to the user during a `filterCells()` function call.
- The `metrics()` function will now look for a stashed `cell2sample` data.frame, which speeds up operations for quality control plots.
- Improved handling of sample metadata columns as factors. In particular, levels should be correctly updated using `droplevels` in a `selectSamples()` call. The [bcbioRNASeq][] package has also been updated to work in a similar fashion, where all columns in the sample metadata data.frame are now defined as factors.
- Simplified `bcbioSingleCell` to `seurat` object coercion to stash all of the bcbio metadata, and simply return the basic `seurat` object, rather than trying to also perform normalization and scaling. These steps have instead been added back to the Seurat R Markdown clustering template.
- Updated cell cycle and cell type markers from our master copy on Google Sheets.
- Added a troubleshooting section to the GitHub README, with a note on maximum DLLs.



# bcbioSingleCell 0.0.21 (2017-11-08)

- Updated package imports to match Bioconductor 3.6.
- Initial support for `plotCellTypesPerCluster()`.
- Initial support for `plotMitoVsCoding()`. I broke this code out from `plotMitoRatio()`. We could opt to keep this in `plotMitoRatio` with a `geom = "scatterplot"` argument.
- Initial commit of `.applyFilterCutoffs()` internal function, used to subset the object to contain only cells and genes that have passed quailty control filtering.
- Improved Seurat `FindAllMarkers()` sanitization.
- Made the quality control plots more modular. Now they support multiple geoms, including `boxplot`, `histogram`, `ridgeline`, and `violin` (default). Median labels are applied with the internal `.medianLabels()` function.
- Updated error message for CellRanger directory structure.
- Cellular barcode columns are no longer split with an underscore. Instead, they are kept as a single ACGT string. We're now generating `cellID` to `sampleID` matching with a different method. In the future, we'll stash a cell2sample data.frame inside the object, that makes this operation faster than the current `mclapply()` code.
- Updated annotable support in `loadSingleCell()` and `loadCellRanger()`.
- Added support for handling both gene- and transcript-level counts. The updated release of the bcbio single-cell pipeline now outputs at gene level.
- Initial quality control plot support for seurat objects.
- Added support for return of only filtered cells and genes with the `counts(filterCells = TRUE)` function.
- Added assignment support for `interestingGroups<-`.
- Initial support for sample metadata generation from seurat object.
- Improved internal code for `selectSamples()`.
- Updated `topMarkers()` to match Seurat v2.1 update.
- Improved `bcbioSingleCell` to `seurat` coercion method with `setAs()`.



# bcbioSingleCell 0.0.20 (2017-10-24)

- Upgraded to basejump 0.1.0 and Seurat 2.1.0 dependencies.
- Improved documentation of NAMESPACE imports per function.
- Switched to base grep functions where applicable (`grepl()`, `gsub()`).
- Use GTF in package documentation rather than GFF. Applies primarily to the `loadSingleCell()` import function.
- Restrict class support in S4 methods to `bcbioSingleCell`. Legacy `bcbioSCDataSet` class can be upgraded to `bcbioSingleCell` class using `as(bcb, "bcbioSingleCell")` coercion.
- Use filtered cell output for `metrics()` and quality control functions by default.
- Updated the quality control R Markdown to include `filterCells = FALSE` where applicable.
- Draft support for aggregated technical replicates in quality control functions using `sampleNameAggregate` column in sample metadata. This doesn't change the actual counts values. It only applies to visualization in the quality control plots currently.
- Miscellaneous R Markdown template updates. Primarily improvements to the setup chunk object loading workflow.
- Removed lintr checks from testthat. This is breaking `devtools::test()`.



# bcbioSingleCell 0.0.19 (2017-10-12)

- Renamed main object class from `bcbioSCDataSet` to `bcbioSingleCell`.
- Cell filtering with `filterCells()` will now slot a named logical vector into `metadata(object)[["filteredCells"]]`, which will be used to dynamically subset the slotted internal `SummarizedExperiment` data. Now that we're using this approach, we can return a modified `bcbioSingleCell` object rather than defining a separate `bcbioSCFiltered` class.
- Renamed `loadSingleCellRun()` to `loadSingleCell()`, to match bcbioRNASeq package.
- Now allowing implicit integers in our function code.
- Added support for plotting technical replicates. This is handled by `sampleNameAggregate` in the sample metadata.
- Now using ridgeline plots in place of histograms where applicable.
- Travis CI checks take too long when loading SummarizedExperiment. Hopefully this will be fixed in the 3.6 release later this month.
- New internal dark theme (`darkTheme()`), based on the Seurat theme.
- Initial commit of `plotDot()` function, based on `Seurat::DotPlot()`.
- Added new tSNE plots that allow for consistent cluster labeling.
- Providing legacy support for `bcbioSCDataSet` and `bcbioSCFiltered`, which will be deprecated in a future release.
- Offloaded some internal code to basejump, for improved consistency with bcbioRNASeq package: `internal-projectDir.R`, `internal-readSampleMetadataFile.R`, `internal-sampleDirs.R`. We may want to provide this code as a shared bcbio core package (e.g. bcbioBase) in the future.
- Added internal utility to check for valid marker genes (`.validMarkers()`).
- Improved Ensembl release version support (`ensemblVersion`).



# bcbioSingleCell 0.0.18 (2017-09-17)

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



# bcbioSingleCell 0.0.17 (2017-09-03)

- Renamed package from `bcbioSinglecell` to `bcbioSingleCell`.
- Added [viridis][] color palette support to quality control plots.
- Added cell-cycle marker genes for human and mouse.
- Slotted `organism` in `bcbioSCDataSet` metadata, in addition to `genomeBuild`.
- Updated `bcbioSCFiltered` to `seurat` coercion method to also run `FindVariableGenes()` and `ScaleData()` by default.
- Added cell-cycle regression into [Seurat][] clustering [RMarkdown][] template.
- Improved [pkgdown][] settings and website appearance.



# bcbioSingleCell 0.0.16 (2017-08-25)

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



# bcbioSingleCell 0.0.15 (2017-08-11)

- Renamed functions in `lowerCamelCase` from `snake_case`.
- Draft support for [monocle][], [scater][], and [scone][].
- Package now depends on [SummarizedExperiment][].
- Renamed `loadRun()` to `loadSingleCellRun()` for improved compatibility with [bcbioRNASeq][] package. This helps avoid NAMESPACE collisions between packages.
- Improved support for custom GTF files.
- Split out [Cell Ranger][] import into a separate utility function named `loadCellRanger()`.
- Added plotZerosVsDepth() by @roryk.
- Renamed `filteringCriteria` to `filterParams` in `@metadata` slot.
- Miscellaneous documentation fixes.



# bcbioSingleCell 0.0.14 (2017-07-26)

- Migrated plotting functions to S4 methods.
- Improved sample metadata handling.
- Updated sample name and cellular barcode sanization.
- Initial commit of [Seurat][] utility functions.
- Initial commit of YAML parameters for [RMarkdown][] templates.
- Draft support for [Seurat][] v2 pre-release.
- Initial support for [SureCell][] UMIs.



# bcbioSingleCell 0.0.13 (2017-07-10)

- Integrated [bcbio][] and [Cell Ranger][] workflows into `load_run()`.



# bcbioSingleCell 0.0.12 (2017-06-28)

- Initial support for S4 object creation using [SummarizedExperiment][].



# bcbioSingleCell 0.0.11 (2017-06-24)

- Added [testthat][], [lintr][], and [covr][] support for code coverage.



# bcbioSingleCell 0.0.10 (2017-06-15)

- Updated [RMarkdown][] templates.



# bcbioSingleCell 0.0.9 (2017-06-12)

- Compatibility update for basejump S4 NAMESPACE changes.
- Improved plots for quality control and filtering.



# bcbioSingleCell 0.0.8 (2017-05-28)

- Renamed mitochondrial plot functions.
- Changed presentation of mitochondrial abundance as ratio instead of percentage.
- Simplified NAMESPACE by offloading dependencies to basejump.



# bcbioSingleCell 0.0.7 (2017-05-13)

- Initial support for loading of 10x Genomics [Cell Ranger][] output.
- Add detection of droplet method based on the metadata file.
- Draft support for import of [inDrop][] i5 index barcode counts.



# bcbioSingleCell 0.0.6 (2017-05-10)

- Modified `load_run()` function for improved consistency with [bcbioRNASeq][] package.
- Initial incorporation of scater package into workflow.
- Initial commit of RMarkdown template skeletons.



# bcbioSingleCell 0.0.5 (2017-05-08)

- Improvements to run object based on code in [bcbioRNASeq][] package.
- Better handling in plots for datasets with many samples.
- Labeled quality control cutoffs on the plots.



# bcbioSingleCell 0.0.4 (2017-04-26)

- Prepared package for [tidyverse][] dependency updates.
- Updated HPC server detection.



# bcbioSingleCell 0.0.3 (2017-04-20)

 - NAMESPACE consolidation and function cleanup.



# bcbioSingleCell 0.0.2 (2017-04-12)

- Added quality control plotting functions.
- Converted all functions to standard evaluation.



# bcbioSingleCell 0.0.1 (2017-03-01)

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
