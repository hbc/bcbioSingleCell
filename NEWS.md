# bcbioSingleCell 0.5.0 (2022-03-11)

## Major changes

- Split out basejump imports into component packages.

## Minor changes

- Removed unused `formalsList` from NAMESPACE, which is removed in upcoming
  basejump package update.
- Removed `barcodeRanksPerSample` from NAMESPACE.

# bcbioSingleCell 0.4.16 (2021-03-12)

## Minor changes

- Updated basejump dependencies and removed unnecessary stringr import.

# bcbioSingleCell 0.4.15 (2021-02-22)

## Minor changes

- Maintenance release, providing support for basejump v0.14 release series.

# bcbioSingleCell 0.4.14 (2020-10-08)

## Minor changes

- Maintenance release, updating dependency version requirements.

# bcbioSingleCell 0.4.13 (2020-07-24)

## Minor changes

- Maintenance release, updating minimum R dependency to 4.0.

# bcbioSingleCell 0.4.12 (2020-05-12)

## Minor changes

- Bug fix for integer handling in `.importReads` following switch to
  `vroom::vroom` from `data.table::fread`.

# bcbioSingleCell 0.4.11 (2020-05-11)

## Minor changes

- `updateObject`: Bug fix release for breaking changes introduced by `colData<-`
  assignment on `SingleCellExperiment` / `SummarizedExperiment` expecting
  `check = FALSE` support for downstream `updateObject` call. Added `...`
  passthrough that prevents this error in method support for `bcbioSingleCell`
  class. This should now be fully compliant with R 4.0 / Bioconductor 3.11
  release.
- `updateObject`: Added `verbose` argument support, defaulting now to `FALSE`.

# bcbioSingleCell 0.4.10 (2020-02-24)

## Minor changes

- `bcbioSingleCell`: Removed now defunct `spikeNames` argument. Refer to
  `makeSummarizedExperiment`, `makeSingleCellExperiment` documentation and
  release notes in basejump package for details.

# bcbioSingleCell 0.4.9 (2020-02-20)

## Minor changes

- Tightened up dependencies in DESCRIPTION in preparation for bioconda release
  update, which breaks with basejump 0.12.0 conda dependency otherwise.

# bcbioSingleCell 0.4.8 (2020-02-19)

## Minor changes

- Changed license from MIT to GPL-3, matching other bcbio and Acid Genomics
  R packages.

# bcbioSingleCell 0.4.7 (2020-01-20)

## Minor changes

- Updated messages to use cli package approach. The same approach is now used in
  the bcbioRNASeq package as well.
- Updated basejump dependencies to support rename from bioverbs to AcidGenerics.

# bcbioSingleCell 0.4.6 (2019-10-30)

## Minor changes

- Updated S4 validity checks and unit tests to support Bioconductor 3.10.
  Internally, we have to adjust checks against `DataFrame` class to also look
  for `DFrame` class, which has changed in latest S4Vectors update. See related
  update in bcbioRNASeq package.
- Full backward compatiblity with Bioconductor 3.9 is currently maintained.
- Updated basejump dependency versions.

# bcbioSingleCell 0.4.5 (2019-09-17)

## Minor changes

- Switched QC template to use `bcb_file` as main file input, similar to
  bcbioRNASeq QC template.
- Updated lintr config.

# bcbioSingleCell 0.4.4 (2019-09-09)

## Minor changes

- `bcbioSingleCell`: Internal generator now calls `importSampleData` using
  `pipeline = "bcbio"` argument, so we don't run into breaking changes when
  handling user metadata in a future basejump update.
- Updated basejump dependency versions.

# bcbioSingleCell 0.4.3 (2019-08-27)

## Minor changes

- NAMESPACE and basejump dependency updates.

# bcbioSingleCell 0.4.2 (2019-08-22)

## Major changes

- `bcbioSingleCell`: Simplifed internal `colData` handling, calling
  `calculateMetrics` after creating `SingleCellExperiment` with
  `makeSingleCellExperiment`. This improves consistency with Chromium package
  and enables better preparation of `rowRanges` containing spike-ins and
  transgenes.

## Minor changes

- Moved `metricsCols` global to basejump, so we can share with Chromium.
- Improved generator consistency with other basejump and bcbio R packages.

# bcbioSingleCell 0.4.1 (2019-08-20)

## Minor changes

- `plotReadsPerCell`: Reworked internal code to use `DataFrame` primarily
  instead of `tbl_df`, so we can remove dplyr dependency.
- Updated basejump and bcbioBase dependency versions.
- Removed some additional unused dependencies, to lighten the package.

# bcbioSingleCell 0.4.0 (2019-08-12)

This update is intended to simplify the package and reduce the amount of
plotting code that is defined to `bcbioSingleCell` class specifically, instead
of `SingleCellExperiment`.

## Major changes

- Offloaded plotting functions to AcidPlots, so we can use some of this code in
  a shared manner with objects imported from 10X Genomics Chromium.

## Minor changes

- Improved documentation consistency by using shared roxygen from AcidRoxygen
  package.
- Updated basejump and bcbio R package dependencies.

# bcbioSingleCell 0.3.18 (2019-07-29)

## Major changes

- Updated R dependency to 3.6. This now requires Bioconductor 3.9+.
- `updateObject`: Simplified method support, dropping `rowRanges` argument.
  No longer attempting to coerce assays back to a simple list, which doesn't
  work well on Bioconductor 3.10 Devel currently.

## Minor changes

- Updated basejump dependencies.
- Improved code coverage over 90%.

# bcbioSingleCell 0.3.17 (2019-07-24)

## Minor changes

- Improved internal S4 function naming consistency.
- Updated basejump dependencies.

# bcbioSingleCell 0.3.16 (2019-07-17)

## Minor changes

- Updated basejump dependencies.
- Improved Travis CI docker configuration.

# bcbioSingleCell 0.3.15 (2019-05-29)

## Minor changes

- Updated `barcodeRanksPerSample` and `plotBarcodeRanks` to support DropletUtils
  update in Bioconductor 3.9.
- Improved Travis and AppVeyor CI configuration.

# bcbioSingleCell 0.3.14 (2019-04-25)

- S4 generic reexport documentation fixes.

# bcbioSingleCell 0.3.13 (2019-04-23)

## Minor changes

- Now importing ggplot2 code from [AcidPlots][] package.
- `plotCellCounts`: Improved default appearance, removing black border.
- Updated basejump and bcbio R package dependencies.

# bcbioSingleCell 0.3.12 (2019-04-18)

## Minor changes

- Deprecated `prepareSingleCellTemplate`, in favor of `prepareTemplate`.
- Switched Travis CI configuration to use `singlecell` Docker image.

# bcbioSingleCell 0.3.11 (2019-04-11)

## Major changes

- Restricted S4 methods back to `bcbioSingleCell` from `SingleCellExperiment`
  where applicable. Generally, this applies to the supported QC plotting
  functions. Previously, the methods were defined against `SingleCellExperiment`
  to also support `CellRanger` S4 class, which has since been moved to the
  [Chromium][] package.

## Minor changes

- Improved code coverage by adding additional unit tests.
- Updated internal reexport method for S4 generics, making `reexports.Rd` file
  no longer necessary.
- Renamed R Markdown QC template from `quality_control` to `quality-control`,
  using kebab case instead of snake case.

# bcbioSingleCell 0.3.10 (2019-04-07)

## Major changes

- `updateObject` is now fully backward compatible for all objects created by
  the package, including versions prior to v0.1.
- Tightened up S4 validity class checks for `bcbioSingleCell`. Similar to
  `bcbioRNASeq` class checks, now requiring these elements in `metadata`:
  `bcbioCommandsLog`, `bcbioLog`, `dataVersions`, `gffFile`, `lanes`,
  `programVersions`, `projectDir`, `runDate`, `wd`, and `yaml`. Previously,
  these weren't required because we were sharing internal validity code checks
  with the `CellRanger` S4 class. `CellRanger` has since been moved to the
  new Chromium package, so it's now appropriate to tighten up bcbioSingleCell.
- Extract method `[` on `bcbioSingleCell` class objects no longer alters the
  version of the object. It only sets `subset = TRUE`, similar to bcbioRNASeq.

## Minor changes

- Bug fix release for formal rename in `emptyRanges` from `mcolsNames` to
  `mcolnames`. This has been changed in [basejump][] 0.10.3.

# bcbioSingleCell 0.3.9 (2019-04-01)

## Minor changes

- Migrated development code to [Acid Genomics][].
- Updated [basejump][] dependency to v0.10 release series.

# bcbioSingleCell 0.3.8 (2019-03-18)

## Minor changes

- Updated basejump and bcbioBase dependencies.
- Updated documentation.
- Added complete release information in NEWS.

# bcbioSingleCell 0.3.7 (2019-02-12)

## Minor changes

- `bcbioSingleCell` generator: improved internal handling of `mapCellsToSamples`
  call. This reworked join step checks the sampleID mapping more carefully.
- Removed `plotGene` from deprecations. This function has been renamed to
  `plotCounts` instead in the basejump package.

# bcbioSingleCell 0.3.6 (2019-01-13)

## Minor changes

- Updated basejump and bcbioBase dependencies, to pass build checks.
- Split out imports into a separate R file.

# bcbioSingleCell 0.3.5 (2018-12-22)

This release reworks the internal assert checks substantially.

## Major changes

- First release that has fully migrated to using goalie package instead of
  assertive package internally for checks.
- Reworked `bcbioSingleCell` S4 validity checks.
- Reworked `bcbioSingleCell` generator function to use cleaner assert checks.
- Migrated to inheriting S4 generics from bioverbs, where applicable.

## Minor changes

- Reworked documentation to use sentence case instead of title case.

# bcbioSingleCell 0.3.4 (2018-12-01)

## Major changes

- `metrics`: Removed code and consolidated `SingleCellExperiment` method support
  in basejump package.

## Minor changes

- Documentation improvements.
- `plotQC`: Improved validity checks.

# bcbioSingleCell 0.3.3 (2018-11-25)

## Minor changes

- Updated basejump and bcbioBase dependencies, to tighten up build checks.
- Documentation and CI configuration fixes.

# bcbioSingleCell 0.3.2 (2018-11-19)

## Major changes

- Migrating `CellRanger` S4 class to separate Chromium R package.
- Migrated additional S4 generics to basejump package (bioverbs).
- Exporting `calculateMetrics` as a user-accessible function.
- Extract (`[`) now returns `bcbioSingleCell` object again.
- `filterCells`: Reworking to return `SingleCellExperiment` instead of
  `bcbioSingleCell` object.
- Improved plotting functions to inherit params from basejump.

## Minor changes

- Removed `barcodePattern` global variable.

# bcbioSingleCell 0.3.1 (2018-11-14)

## Major changes

- Updated NAMESPACE to inherit from basejump instead of any of the basejump
  subpackages, which are subject to rework during this release series.
- Renamed example data from `cellranger_small` to `cellranger`, and
  `indrops_small` to simply `indrops`.

## Minor changes

- Transitioning to goalie for assert checks in place of assertive package.
- Miscellaneous documentation improvements.
- Updated working examples for plotting functions, using `indrops` example data.
- Reworked `bcbioSingleCell` `show` (`print`) method.

# bcbioSingleCell 0.3.0 (2018-11-06)

Fork of `hbc/bcbioSingleCell` code base, for additional development and
improvements prior to Bioconductor submission. The v0.2 release series is
being maintained on the `hbc` organization page for stability.

The v0.3 release series serves to offload some of the code base to the
[basejump][] package, freeing [bcbioSingleCell][] up to be lighterweight and
easier to maintain long-term.

Note that v0.2 release series will remain pinned to basejump v0.7 series.

CellRanger code will be split out into a separate Chromium R package in a future
update, so some of these `SingleCellExperiment`-specific methods need to be
consolidated into the basejump package.

## Major changes

- `aggregateReplicates`, `cell2sample`, `combine`, `plotZerosVsDepth`,
  `sampleData`, `selectSamples`, `subsetPerSample`, `topBarcodes`: Offloaded to
  basejump as `SingleCellExperiment` S4 methods.
- Removed generics offloaded to basejump.
- `barcodeRanksPerSample`: Reworked internal code, using improved
  `countsPerSample` and `ranks` calculations.
- Updated extract method to work on `SingleCellExperiment`, so it can also be
  inherited for `CellRanger` S4 class.
- Reorganized `metricsPerSample` code approach.

## Minor changes

- Consolidated S4 validity check code into `AllClasses.R` file.
- Updated Travis and AppVeyor CI configuration.
- Split out some internal code into `barcodes-internal.R`. This includes
  `.nCount` and `.rawMetrics`.
- Improved documentation for `bcbioSingleCell` generator function.
- Reorganized internal CellRanger processing functions.
- Improved documentation for `CellRanger` generator function.
- Reworked plotting function code to integrate with basejump package.
- Split out plotting code methods into internal functions
  (e.g. `plotMitoVsCoding.SingleCellExperiment`).

# bcbioSingleCell 0.2.1 (2018-08-19)

## Major changes

- Moved clustering and marker templates to [pointillism][] package, since
  all clustering code has been moved there. The quality control template using
  the `filterCells` function is still provided here.

## Minor changes

- `prepareSingleCellTemplate`: Improved internal code to use updated
  `prepareTemplate`, which has been migrated to [basejump][] package.
- Added back the `cellranger_small` example dataset.

# bcbioSingleCell 0.2.0 (2018-08-09)

This is a significant maintenance release designed to simplify the package for
long-term stability. Here we are removing code support for [Seurat][] and
visualization of cell clustering. This functionality has been moved to the new
[pointillism][] package. We have streamlined [bcbioSingleCell][] to simply
import [bcbio][] single-cell data using the `bcbioSingleCell` constructor or
[Cell Ranger][] data with the `readCellRanger` function. The package continues
to provide quality control visualization (see `plotQC`) and low-quality cell
removal with the `filterCells` function.

## Major changes

- Moved clustering and `seurat` object support to [pointillism][]:
  `cellCountsPerCluster`, `cellTypesPerCluster`,
  `clusterCellCountsPerSample`, `fetchData` functions, `plotFeature`
  functions, `plotGene` functions, `plotMarker` functions, `plotPCA`,
  `plotTSNE`, `plotUMAP`.
- Now requiring [bcbioBase][] v0.4+ and [basejump][] v0.7+.

## Defunct functions

- Removed from exports (for simplicifcation): `plotClusters`,
  `plotFeatures`, `quantileHeatmap`, `readMarkers`, `readMarkersFile`.

## Minor changes

- No longer exporting `mapCellsToSamples`. This is now only called internally.
- Updated documentation to use [roxygen2][] v6.1.

# bcbioSingleCell 0.1.18 (2018-07-31)

## Major changes

- Now recommending [BiocManager][] instead of [BiocInstaller][] to install.
- `filterCells` now supports `zinbwave = TRUE` for automatic
  zero-inflation weights calculation, using the [zinbwave][] package.
- Now requiring in `bcbioSingleCell` object validity checks that all
  `sampleData` columns are defined in `colData`, and that the factor levels
  also match. If you run into an error because of this change, post an issue
  on [GitHub][] and we'll sort it out. The validity check method now also checks
  to make sure that all `metricsCols` are defined in `colData`.
- `aggregateReplicates` now keeps track of interesting groups metadata
  better and returns a `SingleCellExperiment` containing this information.
- `seurat` to `SingleCellExperiment` now keeps track of `rowRanges`, if they
  have been stashed inside the `seurat` object.
- Minimal `bcbioSingleCell` (`indrops_small`) and `SingleCellExperiment`
  (`cellranger_small`) objects now contain dimensionality reduction information
  calculated with [Seurat][]. See the scripts in `data-raw/` for how this was
  performed.
- `metrics`: `SingleCellExperiment` method now just returns `colData`
  with `interestingGroups` column defined internally with
  `uniteInterestingGroups`. Previously we updated the `sampleData` to long
  format here dynamically, but this is no longer needed now that `colData`
  always must contain sample data in long format.
- `sampleData<-` assignment method now dynamically updates the sample metadata
  in `colData` to match.
- Simplified `SingleCellExperiment` method extension for `seurat` objects.
- Added support for global [ggplot2][] discrete color palettes using
  `bcbio.discrete.color` and `bcbio.discrete.fill` with `options`.

# Minor changes

- Marker functions are no longer exported as generics, since they only
  currently apply to [Seurat][] and may be deprecated in the future.
- Improved documentation, specifying the accepted classes for each argument.
  This matches the conventions now used in [bcbioBase][] and [bcbioRNASeq][].
- Now using [roxygen2][] v6.1 for documentation.
- `metrics` `matrix` method now uses `rowRanges` instead of `rowData`.
- `mapCellsToSamples` is no longer exported and has been switched to an
  internal function.
- Improved `filterCells` return to include the counts per sample, which
  makes visual inspection and confirmation with `plotCellCounts` easier.
- Added assert checks to plotting functions to stop if
  `interestingGroups = NULL`. Adding support for NULL handling may be a good
  idea in a future release but isn't a priority at the moment.
- Improved factor level handling in `plotCellCounts` and `plotZerosVsDepth`
  to avoid warnings if `colData` and `sampleData` levels don't match
  exactly. This can happen if the sample order is different in `sampleData`
  than in `colData`.

# bcbioSingleCell 0.1.17 (2018-07-21)

## Major changes

- `readCellRanger` now supports HDF5 data (.h5 files) or MatrixMarket Exchange
  Format explicitly, using the `format` argument. The function also now supports
  user-requested import of raw whitelisted cellular barcodes, using the
  `filtered` argument. By default, filtered counts from [Cell Ranger][] are
  imported.
- Removed support for [zingeR][] in `diffExp`. We will consider adding this
  functionality back in a future update, once the package is available on
  [Bioconductor][]. For now, we're recommending usage of [zinbwave][].
- `filterCells` now supports `nCells` argument for hard cutoffs, applied after
  other QC filtering cutoffs using `nUMI`.
- Internal data has been renamed to snake_case format, for improved consistency.
  Note that `cell_cycle_markers` is used in place of `cellCycleMarkers`, and
  `cell_type_markers` in place of `cellTypeMarkers`. There's no easy way to
  link these files as aliases, so this is a necessary breaking change.

## Minor changes

- Improved [R Markdown][] template defaults for quality control and clustering.
- Now importing functions related to differential expression into the package
  NAMESPACE, including [DESeq2][], [edgeR][], and [zinbwave][].
- Parallelization is no longer used to load samples, improving handling in an
  HPC cluster environment. We will consider using [BiocParallel][] in a future
  update.
- Now requiring [ggplot2][] v3.0. Switched to using `aes` instead of
  `aes_string` internally, using new [tidyeval][] syntax.
- Removed additional single-cell toolkit packages in Suggests, including
  [scater][] and [scran][]. This helps speed up build checks on [Travis CI][].
- Reexported functions have been moved to the [bcbioBase][] package.
- Improved internal consistency between the `bcbioSingleCell` and
  `readCellRanger` data import functions.

## Additional notes

- Package still works with [Bioconductor][] 3.6 / [R][] 3.4, but we will be
  migrating to [Bioconductor][] 3.7 / [R][] 3.5 when [conda][] r-base is updated
  to 3.5.1.

# bcbioSingleCell 0.1.16 (2018-06-20)

## New functions

- `cellCountsPerCluster` and `clusterCellCountsPerSample`.

## R Markdown templates

- Now recommending 1000 UMIs by default, and a maximum of 10% mitochondrial
  transcripts.
- [Seurat][] clustering now calculates multiple resolutions, suggesting 0.4,
  0.8, and 1.2 by default. Dimensional reduction plots have been updated to
  support looping of these multiple resolutions.

## Minor changes

- `as(object, "SingleCellExperiment")` coercion now slots stashed metadata
  into `SingleCellExperiment`, if defined.
- Simplified the internal code for `sampleData`. `seurat` objects now share
  the same code as `SingleCellExperiment`, and return `NULL` if the sample
  data is not defined. The `metrics` function continues to slot empty
  metadata in `sampleID`, `sampleName`, and `interestingGroups` if not
  defined.
- `SingleCellExperiment` is used as default method over `seurat` where
  applicable: `fetchData` family, `plotDimensionalReduction` family,
  `plotMarker` family, `plotFeature` family, `plotGene` family,
  and `plotCellTypesPerCluster`. The internal code hasn't changed, it just
  is defined primarily for `SingleCellExperiment`.
- `dimRed` argument has been renamed to `reduction`, where applicable.
- `topBarcodes` can now return either a `data.frame` or `list`, containing
  the top barcodes grouped by sample.

# bcbioSingleCell 0.1.15 (2018-06-13)

## Minor changes

- Updated internal code to use `text` as primary argument in `markdownHeader`
  calls.
- Updated example datasets and unit tests to match.
- Working on using pbmc4k as example dataset for vignette.

# bcbioSingleCell 0.1.14 (2018-06-05)

## Major changes

- Multiplexed [Cell Ranger][] barcodes are now reformatted to be more readable
  and compatible with the `mapCellsToSamples` utility function. This change
  helps simplify the internal code for `cell2sample`. For example in the
  pbmc4k dataset, the barcode IDs are sanitized from `TTTGGTTTCGCTAGCG-1` to
  "`pbmc4k_1_TTTGGTTTCGCTAGCG`". The "`1`" here denotes 1 sample in the matrix,
  which is how [Cell Ranger][] denotes multiplexed samples in a single counts
  matrix. Note that the "`-`" character is illegal in names, so we consistently
  sanitize barcodes to contain "`_`" instead. See `help("make.names")` for
  more information on syntactically valid names in [R][].
- `readCellRanger` no longer requires reference data defined by `refdataDir`,
  although this is still recommended.

## Minor changes

- Resaved example datasets.
- Switched the example `cellranger_small` and `seurat_small` datasets to the
  publicly available pbmc4k dataset from [10X Genomics][]. Here we've subset the
  top 500 cells and genes by abundance. We'll use either the pbmc4k or pbmc8k
  dataset for the vignette in a future update.
- `bcbioSingleCell` and `readCellRanger` functions now consistently default
  to not requiring `sampleMetadataFile`, which is now `NULL` by default. For
  `bcbioSingleCell`, if a custom sample metadata file is not provided, the
  function reads from the bcbio YAML metadata. For `readCellRanger`, the
  function uses `minimalSampleData` internally to return minimal metadata,
  containing `sampleName` and `description` columns.
- Broke out internal `.sampleDirs` function to `bcbioSingleCell` and
  `readCellRanger` functions.
- Improved `plotMarker` documentation examples to use mitochondrial genes.
- Using `seurat_small` in place of `Seurat::pbmc_small` in working examples.

## Internal code changes

- Now consistently using [dplyr][] for piped `data.frame` as much as possible,
  where applicable. Code is being updated to use [tidyeval][].

# bcbioSingleCell 0.1.13 (2018-05-29)

## Minor changes

- t-SNE and UMAP plot improvements.

# bcbioSingleCell 0.1.12 (2018-05-25)

## Major changes

- No longer using automatic camel case sanitization for `metrics` or
  `fetchData` return column names.
- Improved [R Markdown][] clustering and marker templates to optionally support
  UMAP and dark mode in the YAML parameters.

## Minor changes

- Using original [Seurat][] mapping names for data: tSNE_1, tSNE_2, PC1, PC2,
  UMAP1, UMAP2.
- Ensure transcript-level counts always have `stripTranscriptVersions` command
  applied, to remove the [Ensembl][] transcript versions if present.
- No longer using labels (e.g. A, B, C, D) on `ggplot` grid return.
- Note that "phase" has been renamed to "Phase" in the [R Markdown][] clustering
  for cell-cycle regression PCA.

# bcbioSingleCell 0.1.11 (2018-05-23)

## New functions

- UMAP is now supported. This functionality is provided in: `plotUMAP`,
  `plotMarkerUMAP`, and `plotFeatureUMAP`. Corresponding fetch functions,
  `fetchUMAPData` and `fetchUMAPExpressionData`, have also been added.
- `plotGene`: Added `seurat` method support. If advanced customization of
  the plot is needed, use `plotDot` or `plotViolin` instead, or refer to the
  [Seurat][] documentation for alternates.

## Major changes

- Dimensional reduction and marker plots no longer use dark mode by default.
  The default color palette support for marker plots has been improved to
  consistently use [viridis][].
- `diffExp`: improved internal code to work directly on
  `SingleCellExperiment`, removing the need to pass `design` and `group`
  parameters internally. Also added unit testing against [zinbwave][],
  [zingeR][], and [edgeR][] support. [DESeq2][] is supported but runs slowly.
- Reworked `plotFeature` and `plotMarker` family of functions. Improved the
  color palette support when `dark = FALSE`, now using a flipped [viridis][]
  plasma color palette.
- `aggregateReplicates` function has been reworked to return a
  `SingleCellExperiment` object instead of `bcbioSingleCell`. The v0.2.4
  update of [bcbioRNASeq][] behaves similarly with this generic.

## Minor changes

- Reworked the internal handling of some `seurat` `SingleCellExperiment` method
  support, using `as(x, "SingleCellExperiment")` internally, which uses the
  new `Seurat::Convert` function.
- Made some previously deprecated functions now defunct: `plotClusters`,
  `plotTSNEExpressionData`, `loadSingleCellRun`, `darkTheme`,
  `pcCutoff`, `quantileHeatmap`, `plotKnownMarkers`, `readMarkers`,
  `readMarkersFile`.
- Made `plotFeatures`, `plotMarker`, and `plotMarkers` functions defunct.
- `plotPCElbow` now returns a plot grid.
- `sanitizeMarkers`: improved internal code for supported bcbio stashed
  metadata, including `rowRanges`.

# bcbioSingleCell 0.1.10 (2018-05-19)

## Minor changes

- `plotCellTypesPerCluster` is using `dark = TRUE` by default again.
- Fixed `cell2sample` handling for multiplexed [Cell Ranger][] data loaded up
  with `readCellRanger`. Need to use stashed `cell2sample` factor saved in
  `metadata`, rather than attempting to calculate on the fly with
  `mapCellsToSamples`.
- Updated [Travis CI][] build checks to include bioc-release on [macOS][].

# bcbioSingleCell 0.1.9 (2018-05-18)

## Major changes

- No longer attempting to sanitize the rownames for `seurat` objects in coercion
  method. This helps maintain the gene symbol appearance in plotting functions
  for genes with hyphens in the names.

## Minor changes

- Using `BiocParallel::SerialParam` internally for zinbwave in `diffExp`.
- Simplified `cell2sample` internal code to always use `mapCellsToSamples`
  instead of attempting to use a stashed `vector` inside `metadata` for
  `SingleCellExperiment` method.
- Removed internal `.applyFilterCutoffs`, which is no longer necessary since
  this functionality is supported in the S4 subset method.
- Simplified assert checks inside `fetchGene` functions.
- `plotCellTypesPerCluster`: revert back to `dark = TRUE` by default.
- Consolidated `plotMarker` and `plotFeature` functions in the documentation.
- `sanitizeMarkers`: Improved gene identifier matching.
- `topMarkers` now defaults to `coding = FALSE` by default, since not all
  datasets will contain biotype information.

# bcbioSingleCell 0.1.8 (2018-05-16)

## Minor changes

- Initial `updateObject` method support for `bcbioSingleCell` class.
- Relaxed `validObject` validity check to not require sample-level metadata in
  `colData` yet.

# bcbioSingleCell 0.1.7 (2018-05-15)

## New functions

- Added support for Uniform Manifold Approximation and Projection (UMAP) with
  the `plotUMAP` and `fetchUMAPData` functions. These work similarly to the
  other `plotDimensionalReduction` and `fetchData` functions.

## Major changes

- Now adding sample-level metadata into `colData` slot, for better downstream
  compatibility with other packages that work with `SingleCellExperiment`
  container class. Unique per-sample rows are still saved internally in the
  `sampleData` slot.
- Now recommending ECDF as the default geom for quality control plots, where
  applicable.
- `filterCells` now supports `minUMIs = c("knee", "inflection")` for automatic
  filtering based on the cellular barcode ranks. Internally this is handled
  by `DropletUtils::barcodeRanks`.

## Minor changes

- Attempting to re-enable `libgsl-dev` installation for [zinbwave][] on
  [Travis CI][].
- Suggesting [BiocParallel][] for [zinbwave][] call in `diffExp`.
- Now importing [Seurat][] functions into NAMESPACE.
- Consolidated `fetchData` functions in the documentation.
- Consolidated `plotDimensionalReduction` functions in the documentation.
- Updated `aggregateReplicates` internal code. This function again only
  supports aggregation of `bcbioSingleCell` objects that have been filtered
  using the `filterCells` function.
- Now using `Seurat::Convert` internally to coerce `seurat` class object
  to `SingleCellExperiment`, using `as(seurat, "SingleCellExperiment")`. This
  utility function was added to [Seurat][] v2.3.1.

# Internal changes

- Tweaked `metrics` `SingleCellExperiment` method code to always merge
  `colData` and `sampleData`.
- Updated `readCellRanger` internal code to match `bcbioSingleCell`
  constructor, specifically handling sample-level metadata in `colData`.

# bcbioSingleCell 0.1.6 (2018-05-09)

## Minor changes

- Updated default QC [R Markdown][] template.
- Added trendline option to QC scatter plot functions.
- Simplified internal handling of interestingGroups in `plotQC`.
- Using `readYAMLSampleData` internally instead of defunct
  `sampleYAMLMetadata`.
- Added `sampleNames` method support for seurat.
- Now importing rmarkdown, sessioninfo, tidyverse for R Markdown reports, rather
  than suggesting. Similar update applied to bcbioRNASeq.
- Suggesting scater and scran.

# bcbioSingleCell 0.1.5 (2018-05-04)

## Major changes

- Overhauled inflection and knee point labeling support in `plotUMIsPerCell`.
  Now uses the `point` argument and always labels per sample. Currently requires
  the `geom = "ecdf"` argument for labeling.
- Updated default quality control template.
- Added `plotBarcodeRanks`.

## Minor changes

- Added barcode rank support for `seurat` class objects.
- QC plots now have titles by default, matching the conventions used in
  bcbioRNASeq.
- Fixed y-axis scale for histogram geom in QC plots.
- Prefiltering of very low quality barcodes with no UMIs or genes is now always
  applied. This helps avoid unwanted downstream errors with zero count barcodes.
- Added boxplot geom support for `plotReadsPerCell()`.
- `plotQC` geom argument is now more consistent across the paneled plots.
- Fixed facet wrapping for aggregate samples in the QC plots.
- Added `interestingGroups` support to `plotZerosVsDepth`, matching the other
  QC functions.

# bcbioSingleCell 0.1.4 (2018-04-30)

- Updated `sampleData` S4 methods to match update in bcbioBase. Now supports
  `clean` argument, which returns non-blacklisted factor columns only. See
  `bcbioBase::metadataBlacklist` for the blacklist.
- Improved axis scale appearance on dimensionality reduction plots using
  `scales::pretty_breaks` internally.
- Added `grid` argument to plots, where applicable.
- Renamed example dataset from `bcb_small` to `indrops_small`.
- Removed unnecessary method support for `interestingGroups` and `metadata`.
  These extend from `SummarizedExperiment` correctly now.
- Fixed x-axis label centering for `plotCellCounts` and `plotReadsPerCell`.
- Simplified seurat method support for `SingleCellExperiment`-like methods,
  where applicable. This includes `rowData`, `gene2symbol`, and
  `interestingGroups`.
- Improved dark mode color support for `plotDot`, `plotFeatureTSNE`,
  `plotMarker`.
- Updated sample metadata example in README.

# bcbioSingleCell 0.1.3 (2018-04-25)

## Minor changes

- Improved summary statistics output during `filterCells` call.
- Miscellaneous documentation improvements, most notably to `bcbioSingleCell`
  constructor function.
- `plotViolin` now uses a color border by default.
- Improved `cell2sample` mapping internally for `readCellRanger`.
- Improved [Bioconductor][] 3.7 installation instructions.

# bcbioSingleCell 0.1.2 (2018-04-24)

## Major changes

- Now using `bcbioSingleCell` instead of `loadSingleCell` as the main
  constructor function to create a `bcbioSingleCell` object. `loadSingleCell`
  is deprecated and still works, but will warn the user.
- Renamed `loadCellRanger` to `readCellRanger` for better name consistency.
- Quality control function color palettes now default to [ggplot2][] colors
  instead of using [viridis][] palettes. This is defined using
  `scale_color_hue` instead of `scale_color_viridis` for example. The
  [viridis][] color palette is still used by default for marker expression
  plots.
- Use "`aggregate`" instead of "`sampleNameAggregate`" to define
  aggregate/grouped samples in metadata.

## Minor changes

- Reexporting relevant [ggplot2][] and [viridis][] color palettes.
- Renamed references to "inDrops" from "inDrop", where applicable.
- Consistently using `sym` in place of `.data` internally for tidy code.
- Updated `seurat` blacklist for `sampleData` generic.
- `plotCellTypesPerCluster` and `plotMarkerTSNE` now use an automatic color
  palette by default, which enables for dynamic color palette support when
  `dark = TRUE`. Internally this is handled with the `theme_midnight` and
  `theme_paperwhite` [ggplot2][] themes.
- Updated installation instructions to support [Bioconductor][] 3.7.

## Internal changes

- `metrics` method support now defaults to `matrix` and works similarly
  for `dgCMatrix` sparse matrices. This is used in place of `calculateMetrics`
  to generate the per cell quality control metrics.
- Reworked and improved `aggregateReplicates` internal code.
- Fixed facet wrapping when `aggregate` is defined in metadata for quality
  control plots.
- Improved internal code for `sanitizeMarkers` to use map the gene annotations
  from `rowRanges` better.

# bcbioSingleCell 0.1.1 (2018-04-16)

## Major changes

- Added support for calculating `barcodeRanks` and `barcodeRanksPerSample`.
- Now exporting `plotMarker` in addition to `plotMarkers`.
- Primary counts matrix slotted into `assay` is named `counts` instead of
  `raw`, for better consistency with `SingleCellExperiment` class. The
  `counts` generic requires that the primary assay slot is named `counts` to
  work correctly. Nothing else here has changed, just the name.
- `loadCellRanger` now returns a `SingleCellExperiment` object instead of a
  `bcbioSingleCell` object.
- Added support for `transgeneNames` and `spikeNames` when loading up a dataset.
- Datasets from a poorly annotated genome can now be loaded up using
  `organism = NULL` during the `loadSingleCell` call.
- Methods now dispatch on `SingleCellExperiment` rather than `bcbioSingleCell`
  where applicable, providing support for `SingleCellExperiment` objects created
  elsewhere.

## Minor changes

- Added support for labeling barcode ranks, such as elbow or inflection point on
  UMI counts per cell plots.
- Added support for plotting metrics using an empirical distribution function
  (ECDF) plot.
- Use `theme_midnight` and `theme_paperwhite` internally for dimensionality
  reduction plots.
- TSNE and PCA plots now use an aspect ratio of 1 by default.
- Improved imports from Matrix, S4Vectors, and methods.
- Switched back to using base `stop`, `warning`, and `message`.
- `inflectionPoint` has been made defunct, in favor of using `barcodeRanks`.
- Moved internal constructors into S4 methods, where applicable.

# bcbioSingleCell 0.1.0 (2018-04-04)

## Major changes

- `bcbioSingleCell` S4 class now extends `SingleCellExperiment` instead of
  `SummarizedExperiment`. This requires definition of `rowRanges` inside the
  object instead of `rowData`. Similar functionality was added to the
  bcbioRNASeq package. Upgrade support will be provided using `updateObject`in
  a future release.
- Added a differential expression utility function named `diffExp`, which uses
  [zingeR][]/[edgeR][] internally to calculate gene expression changes across
  cell groups.

## Minor changes

- Added `plotCumulativeUMIsPerCell` utility. This may be removed in a future
  update in favor of adding this plot into `plotUMIsPerCell` using an ECDF
  plot.
- Now using `readCellTypeMarkers` to load marker `data.frames`, instead of
  `readCellTypeMarkersFile`. This matches the conventions used in the
  [bcbioBase][] package.

# bcbioSingleCell 0.0.32 (2018-03-12)

- Fixes for object subsetting and `loadSingleCell` organism calls

# bcbioSingleCell 0.0.31 (2018-02-21)

- Improved handling of old [Ensembl][] release version for [Cell Ranger][]
  output.
- Switched [R Markdown][] templates to consistent snake case formatting for
  variables and file names.
- `prepareSingleCellTemplate` now uses `_setup.R` instead of `setup.R`.
- Updated unit tests to work with new assert checks.
- Export `mapCellsToSamples` instead of `cell2sample`. `cell2sample` now
  simply acts as an accessor function, returning the internally stored
  cell2sample mappings rather than trying to calculate. `mapCellsToSamples`
  performs the actual mapping from cellular barcodes to sample identifiers.
- Export `metricsPerSample`
- Now using `BiocParallel::bpmapply` to loop across the sparse matrix files
  per sample internally in the `.sparseCountsList` function, which is shared
  between `loadSingleCell` and `loadCellRanger`.
- Updated assert checks.
- Updated `prepareSingleCellTemplate` to explicitly state which files to
  include for each R Markdown template, rather than inheriting from
  `bcbioBase::prepareTemplate`.
- `selectSamples` now fails on a sample mismatch, rather than warning.
- Improved internal code for subset method, adding a `drop = FALSE` argument to
  the cellular barcode matching call.
- Made all file loads in working examples a single line.
- Simplified [R Markdown][] directory paths.
- Made all loads in unit tests a single line.
- Renamed `programs` metadata slot to `programVersions`, to improve consistency
  with [bcbioRNASeq][] package.
- Updated internal sample metadata sanitization to ensure all columns are
  factors.

# bcbioSingleCell 0.0.30 (2018-01-26)

- Fixes for QC plot labeling of bcbioSingleCell object with different per sample
  filtering cutoffs applied.
- Added new `fetchGeneData` function, that wraps the functionality of
  `Seurat::FetchData` for specific genes.
- Improved internal dimensionality reduction plotting code. The main function
  has been renamed to `.plotDR` internally.
- Simplified the code for `fetchTSNEExpressionData` to return a standard
  `data.frame` with the cellular barcodes as rows, instead of the previously
  grouped `tibble` method. Now this function returns aggregate gene marker
  calculations in the `mean`, `median`, and `sum` columns. Since this method has
  way fewer rows than the grouped `tibble`, the [ggplot2][] code for
  `plotMarkerTSNE` now runs faster.
- Explicitly declare `viridis::` for color palettes, where applicable.
- `colorPoints` argument has been renamed to `expression` for
  `plotMarkerTSNE`.
- Dynamic gene symbol to ensgene conversion has been removed from the plotting
  functions, for greater simplicity. Now the `genes` argument simply matches
  against the rownames in the counts matrix of the object.
- Decreased the default `minCumPct` argument for `plotPCElbow` from 0.9 to
  0.8. This is more conservative and will return slightly fewer principal
  components for dimensionality reduction, by default.
- `plotTSNE`, `plotPCA`, and the other dimensionality reduction-related
  plotting functions now default to a smaller point size (0.5) and slight alpha
  transparency (0.8), to make super imposed points more obvious for large
  datasets with many cells.
- Added a working example for `subsetPerSample`.
- Internally switched from `.onLoad` to `.onAttach` method for automatically
  loading required dependency packages.
- `plotPCA` now uses `phase` instead of `Phase` plotting cell cycle regression
  as an interestingGroup (see Seurat clustering template). Previously some of
  the [Seurat][] metadata columns were not consistently sanitized to
  lowerCamelCase (e.g. `Phase`, `res.0.8`, `orig.ident`).
- Suppress package startup messages in [R Markdown][] templates.

# bcbioSingleCell 0.0.29 (2018-01-24)

- Switched to [rlang][] methods for errors, messages, and warnings: `abort`,
  `inform`, and `warn`.
- Updated `filterCells` function to enable per sample filtering cutoffs. This
  works by passing in a named numeric vector, where the names must match the
  internal `sampleID` metadata column (not `sampleName`).
- Improved internal sanitization of metrics available with the `metrics`
  accessor. Now all count columns (e.g. `nUMI`) are consistently integers, and
  all character vector columns are consistently coerced to factors.
- Seurat metadata available through the bcbioSingleCell generics are now
  consistently sanitized in lowerCamelCase. This applies to `orig.ident` and the
  `res.*` metadata columns.
- Explicit integers are now consistently used in all of the function parameter
  arguments.

# bcbioSingleCell 0.0.28 (2018-01-22)

- Manually define functions used to read barcodes and matrices. This improves
  the functionality of the internal `.readSparseCounts` function.
- Improved sample directory matching for [Cell Ranger][] output.
- Fixed sample metadata subsetting based on cell2sample factor levels in the
  internal subset code.
- Added tabbed histogram, violin, and barplots for appropriate quality control
  functions in the [R Markdown][] code.

# bcbioSingleCell 0.0.27 (2018-01-18)

- Migrated core dependency imports from [basejump][] to [bcbioBase][].
- Improved `colnames` and `rownames` handling for internal `.readSparseCounts`
  function.
- Reworked `loadCellRanger` function. `refDataDir` parameter has been renamed
  to `refdataDir`.
- Added `organism` and `genomeBuild` options to `loadSingleCell`, to override
  the metadata set in the bcbio run, if necessary.
- Improved if statement data class checks, where applicable.
- Renamed internal `.sparseCountsTx2Gene` to `.transcriptToGeneLevelCounts`.
- Updated `geomean` bind method in `fetchTSNEExpressionData`.
- Changed `minNovelty` default from 0.8 to 0.75.
- Improved seurat class support for `plotDot`.
- Improved parameter names for `plotKnownMarkersDetected`. Now uses
  `tsneColor`, `violinFill`, and `dotColor`. Also added `pointsAsNumbers`
  parameter.
- Added `subtitle` parameter for `plotMarkerTSNE`.
- Added `tsneColor`, `violinFill`, `dotColor`, and `dark` parameters for
  `plotMarkers`.
- Improved looping method for `plotTopMarkers` so that it renders correctly in
  [R Markdown][] calls.
- Added method support for `plotViolin`.
- Improved internal `cell2sample` handling in subset method code.

# bcbioSingleCell 0.0.26 (2017-12-18)

- Renamed `readMarkersFile` to `readCellTypeMarkersFile`.
- Improved internal handling of multiplexed CellRanger samples. These are count
  matrices with cellular barcodes ending in `-2`, for example.
- Updated `aggregateReplicates` code to work with basejump generic, which uses
  `groupings` instead of `cells` as the grouping parameter.
- Added method support for `detectOrganism`.
- Improved filtering parameter output in `filterCells`.
- Updated `gene2symbol` method support for bcbioSingleCell and seurat objects.
- Fixed working example in `knownMarkersDetected`.
- Reworked internal code for `plotCellTypesPerCluster`.
- Improved internal checks for facet wrapping and color palette parameter
  arguments, where applicable.
- Offloaded base `plotQuantileHeatmap` functionality into basejump, for use in
  [bcbioRNASeq][] package.
- Added unit test data for `loadSingleCell()`.
- Added additional unit tests to improve code coverage.

# bcbioSingleCell 0.0.25 (2017-12-11)

- Prepared a minimal working example dataset.
- Added some initial unit tests.
- Renamed `plotKnownMarkers` to `plotKnownMarkersDetected`.
- Renamed `readMarkers` to `readMarkersFile`.
- Added `gene2symbol` method support.
- Moved `plotDot` generic to basejump package.
- Added internal `.checkFormat` function, which will check for `ensgene` or
  `symbol` input.
- Added internal `.convertGenesToSymbols` utility function, for mapping
  Ensembl gene identifiers to gene symbols.
- Use explicit calls to Seurat functions interally, for clarity.
- Simplified bcbioSingleCell object return in `loadSingleCell` function.
- Added seurat class method support for `annotable` function.
- Added `bcbio<-` assignment method support for seurat class objects.
- `calculateMetrics` function now uses `annotable = TRUE` as default, instead
  of using `missing` method.
- Added method support for `cell2sample` for seurat class objects.
- `cellTypesPerCluster` now uses `min` and `max` default arguments that don't
  remove any rows.
- Added method support for seurat class objects to `counts` function. This
  defaults to returning the raw counts (`normalized = FALSE`), but can also
  return log-normalized (`normalized = TRUE`) and scaled
  (`normalized = "scaled"`) counts.
- Added Ensembl gene identifier mapping support to `fetchTSNEExpressionData`.
- `filterCells` now simply works in a destructive manner. We've removed the
  `drop` parameter. The messages displayed to the user during this function call
  have been improved, and now include more statistics on the step where the
  majority of cells are filtered.
- Initial commit of `gene2symbol` method support for bcbioSingleCell and
  seurat class objects.
- Improved `interestingGroups<-` assignment method support for seurat class
  objects.
- Color palettes now default to viridis instead of inferno palette, where
  applicable.
- Added Ensembl gene identifier support for `plotDot` function.
- `plotFeatureTSNE` now uses a plural `features` parameter instead of
  `feature`, which is consistent with the syntax used in the other functions.
- Added Ensembl gene identifier support to `plotMarkerTSNE`. The `format`
  argument still defaults to "symbol", for consistency with previous behavior.
  However, in the future we recommend that users pass in stable Ensembl gene
  identifiers here if possible, for better reproducibility.
- `plotMarkers` now supports Ensembl gene identifiers.
- `plotPCElbow` now silently returns the sequence of principal components
  (PCs) that we recommend to use for dimensionality reduction.
- We're now using "glm" instead of "gam" for `geom_smooth` plotting, where
  applicable. See the `plotQC` function code.
- Improved the defaults for `plotQuantileHeatmap` to enable faster plotting.
  Now the dendrogram calculations are skipped by default, which take a long time
  for large datasets.
- Removed the draft `plotStressGenes` function for now. Will try to add this
  in a future update.
- Fixed Markdown header handling for `plotTopMarkers` if `headerLevel = NULL`.
- Simplified `selectSamples` code to rely upon output of our bracket-based
  subsetting method. See `subset.R` file for more details.
- Improved metadata update in bracket-based subsetting method.
- `subsetPerSample` function now defaults to saving in the working directory.
- Updated `topBarcodes` function to rank by `nUMI` instead of `nCount` column,
  so it works with data from either bcbioSingleCell or seurat objects.
- Minor tweaks to seurat coercion from bcbioSingleCell object. Here we've
  improved the metadata slotting.
- Initial commit of example data script in the `data-raw/` directory!

# bcbioSingleCell 0.0.24 (2017-11-27)

- Raw cellular barcodes are now slotted in `object@cellularBarcodes` as a
  `data.frame` instead of a per sample `list`. This makes downstream subsetting
  operations on the barcodes simpler.
- Bug fixes for `cell2sample` mapping.
- Switched back to stable CRAN version of roxygen2 (6.0.1) for documentation.
- Renamed `pcCutoff` to `plotPCElbow`. The function now returns a PC
  sequence from 1 to the cutoff (e.g. 1:10) instead of just the final PC cutoff
  value. The R Markdown clustering template has been updated to reflect this
  change.
- Renamed `quantileHeatmap` to `plotQuantileHeatmap`, for consistency with
  other plotting functions.
- Initial support for customized bracket based subsetting, which now acts upon
  the raw cellular barcode counts stashed in the `object@bcbio` slot.
- Moved `darkTheme` to [basejump][] package and reworked as `midnightTheme`,
  with improved colors and axis appearance.
- Added `pointsAsNumbers` parameter to `plotTSNE` and `plotPCA` functions,
  to match the functionality in `plotMarkerTSNE`.
- Overhauled `loadCellRanger` to support multiplexed [Cell Ranger][] matrix
  output. [Cell Ranger][] adds a numeric suffix to the end of multiplexed
  barcodes (e.g. `AAACCTGGTTTACTCT-1` denotes cellular barcode
  `AAACCTGGTTTACTCT` is assigned to sample `1`).
- Improved `cell2sample` mapping in `aggregateReplicates` function, which uses
  the `sampleNameAggregate` column in sample metadata to define the aggregate
  sample pairings. The `summarize` step at line 101 is slow for datasets with
  many samples and should be changed in the future to speed things up.
- Improved internal `cell2sample` code to handle `NULL` stashed mappings
  better.
- Updated TNSE plotting functions to use `midnightTheme` instead of
  `darkTheme`.
- Added user-defined point and label sizes for `plotMarkerTSNE`.
- Fixed typo in `plotMitoRatio` where `maxGenes` cutoff was plotted instead of
  `maxMitoRatio`.
- Added `legend` parameter argument to `plotQC` function. Also improved
  handling of `NULL` return for `plotReadsPerCell`, which can happen with
  [Cell Ranger][] output.
- Updated facet wrapping in `plotZeroesVsDepth` to match the behavior in the
  other plotting functions.
- Initial methods support for custom bracket-based subsetting.

# bcbioSingleCell 0.0.23 (2017-11-22)

- Improved facet wrapping of aggregated samples (`sampleNameAggregate` present
  in sample metadata), but removing code support for wrapping by multiplexed
  FASTQ description.
- Simplified handling of `bcbioSingleCell` objects with `filterCells` applied.
  This information is stored in the `metadata` slot as 3 variables: (1)
  `filterParams`, numeric vector of the parameters used to define the cell
  filtering cutoffs; (2) `filterCells`, character vector of the cellular barcode
  IDs that passed filtering; (3) `filterGenes`, character vector of the
  [Ensembl][] gene identifiers that have passed filtering, as determined by the
  `minCellsPerGene` parameter.
- For `filterCells` return, we're now defaulting to a destructive operation,
  where the columns (cells) and rows (genes) of the object are adjusted to match
  the cells and genes that have passed filtering. Currently this can be adjusted
  with the `drop` argument for testing, but should generally be left as
  `drop = TRUE`.
- We're now slotting a `cell2sample` named factor in the `metadata` slot,
  which makes downstream quality control operations faster. This is generated on
  the fly for previously saved objects that don't have a stashed `cell2sample`.
- Initial commit of `plotQC` utility function, which plots multiple quality
  control functions for easy visualization. This defaults to output as a cowplot
  grid (`return = "grid"`), but can alternatively be set to return
  [R Markdown][] code (`return = "markdown"`).
- `aggregateReplicates` operation has been improved to properly slot raw
  cellular barcodes in `object@bcbio$cellularBarcodes`. The `filterCells` vector
  is adjusted, and `sampleMetadata` factors should be properly releveled.
- The `counts` accessor simply returns the sparse matrix contained in the
  `assay` slot. The `filterCells` argument has been removed.
- Messages have been added to the `filterCells` function, to help the user
  determine at which step the majority of cells are being filtered. We're
  keeping a non-destructive option using `drop = FALSE` for the time being, but
  this will likely be removed for improved simplicity in a future update.
- Updated the internal code for `metrics` to use a simpler join operation on
  the `colData`, `cell2sample` and `sampleMetadata`.
- Updated facet wrap code in quality control plots to not facet multiplexed
  FASTQ descriptions and simply check for `sampleNameAggregate`.
- Improved appearance of `plotReadsPerCell` labels and legends. Additionally,
  `plotReadsPerCell` more efficiently handles the stashed values in the
  `nCount` column of `colData`, for faster plotting that having to rely on
  manipulation of the raw `cellularBarcodes` list stashed in
  `object@bcbio$cellularBarcodes`.
- `sampleMetadata` return is now consistently sanitized for `bcbioSingleCell`
  and `seurat` objects.
- Minor tweaks to quality control template and setup.R files.
- Added `plotFeatureTSNE` utility function. This improves on
  `Seurat::FeaturePlot` and enables the user to overlay the cluster
  identifiers on top of the t-SNE plot. `plotFeatures` is now deprecated in
  favor of this function.
- Improved internal handling of Seurat data in the `.fetchDimDataSeurat`
  function. This now keeps the cell ID as the rowname.
- Allow the user to define the color palette (`color`), as well as `pointSize`
  and `labelSize` for `plotPCA` and `plotTSNE`.
- Factors are now correctly releveled in `cell2sample` return.
- Improved internal code for `fetchTSNEExpressionData`.
- Bug fix for `metrics` accessor not including the cell ID as rownames.
- Clustering template fixes. Now uses `plotFeatureTSNE` to assess quality
  control metrics on t-SNE.

# bcbioSingleCell 0.0.22 (2017-11-17)

- Now internally stashing a `cell2sample` data.frame, which helps speed up
  operations on cellular barcode metrics calculations for quality control plots.
- Improved support for optional `annotable`, `ensemblVersion`, `gtfFile`, and
  `sampleMetadataFile` arguments in `loadSingleCell` function.
- Simplified some of the messages shown during sample loading, in an attempt to
  make them clearer and more informative. Also improved messages shown to the
  user during a `filterCells` function call.
- The `metrics` function will now look for a stashed `cell2sample`
  `data.frame`, which speeds up operations for quality control plots.
- Improved handling of sample metadata columns as factors. In particular, levels
  should be correctly updated using `droplevels` in a `selectSamples` call.
  The [bcbioRNASeq][] package has also been updated to work in a similar
  fashion, where all columns in the sample metadata data.frame are now defined
  as factors.
- Simplified `bcbioSingleCell` to `seurat` object coercion to stash all of the
  bcbio metadata, and simply return the basic `seurat` object, rather than
  trying to also perform normalization and scaling. These steps have instead
  been added back to the Seurat R Markdown clustering template.
- Updated cell cycle and cell type markers from our master copy on
  [Google Sheets][].
- Added a troubleshooting section to the GitHub README, with a note on maximum
  DLLs.

# bcbioSingleCell 0.0.21 (2017-11-08)

- Updated package imports to match [Bioconductor][] 3.6.
- Initial support for `plotCellTypesPerCluster`.
- Initial support for `plotMitoVsCoding`. I broke this code out from
  `plotMitoRatio`. We could opt to keep this in `plotMitoRatio` with a
  `geom = "scatterplot"` argument.
- Initial commit of `.applyFilterCutoffs` internal function, used to subset
  the object to contain only cells and genes that have passed quailty control
  filtering.
- Improved [Seurat][] `FindAllMarkers` sanitization.
- Made the quality control plots more modular. Now they support multiple geoms,
  including `boxplot`, `histogram`, `ridgeline`, and `violin` (default). Median
  labels are applied with the internal `.medianLabels` function.
- Updated error message for [Cell Ranger][] directory structure.
- Cellular barcode columns are no longer split with an underscore. Instead, they
  are kept as a single ACGT string. We're now generating `cellID` to `sampleID`
  matching with a different method. In the future, we'll stash a cell2sample
  data.frame inside the object, that makes this operation faster than the
  current `mclapply` code.
- Updated annotable support in `loadSingleCell` and `loadCellRanger`.
- Added support for handling both gene- and transcript-level counts. The updated
  release of the bcbio single-cell pipeline now outputs at gene level.
- Initial quality control plot support for seurat objects.
- Added support for return of only filtered cells and genes with the
  `counts(filterCells = TRUE)` function.
- Added assignment support for `interestingGroups<-`.
- Initial support for sample metadata generation from seurat object.
- Improved internal code for `selectSamples`.
- Updated `topMarkers` to match Seurat v2.1 update.
- Improved `bcbioSingleCell` to `seurat` coercion method with `setAs`.

# bcbioSingleCell 0.0.20 (2017-10-24)

- Upgraded to [basejump][] 0.1.0 and [Seurat][] 2.1.0 dependencies.
- Improved documentation of NAMESPACE imports per function.
- Switched to base grep functions where applicable (`grepl`, `gsub`).
- Use GTF in package documentation rather than GFF. Applies primarily to the
  `loadSingleCell` import function.
- Restrict class support in S4 methods to `bcbioSingleCell`. Legacy
  `bcbioSCDataSet` class can be upgraded to `bcbioSingleCell` class using
  `as(bcb, "bcbioSingleCell")` coercion.
- Use filtered cell output for `metrics` and quality control functions by
  default.
- Updated the quality control R Markdown to include `filterCells = FALSE` where
  applicable.
- Draft support for aggregated technical replicates in quality control functions
  using `sampleNameAggregate` column in sample metadata. This doesn't change the
  actual counts values. It only applies to visualization in the quality control
  plots currently.
- Miscellaneous [R Markdown][] template updates. Primarily improvements to the
  setup chunk object loading workflow.
- Removed lintr checks from testthat. This is breaking `devtools::test`.

# bcbioSingleCell 0.0.19 (2017-10-12)

- Renamed main object class from `bcbioSCDataSet` to `bcbioSingleCell`.
- Cell filtering with `filterCells` will now slot a named logical vector into
  `metadata(object)[["filteredCells"]]`, which will be used to dynamically
  subset the slotted internal `SummarizedExperiment` data. Now that we're using
  this approach, we can return a modified `bcbioSingleCell` object rather than
  defining a separate `bcbioSCFiltered` class.
- Renamed `loadSingleCellRun` to `loadSingleCell`, to match [bcbioRNASeq][]
  package.
- Now allowing implicit integers in our function code.
- Added support for plotting technical replicates. This is handled by
`sampleNameAggregate` in the sample metadata.
- Now using ridgeline plots in place of histograms where applicable.
- Travis CI checks take too long when loading SummarizedExperiment. Hopefully
  this will be fixed in the 3.6 release later this month.
- New internal dark theme (`darkTheme`), based on the Seurat theme.
- Initial commit of `plotDot` function, based on `Seurat::DotPlot`.
- Added new t-SNE plots that allow for consistent cluster labeling.
- Providing legacy support for `bcbioSCDataSet` and `bcbioSCFiltered`, which
  will be deprecated in a future release.
- Offloaded some internal code to basejump, for improved consistency with
  [bcbioRNASeq][] package: `internal-projectDir.R`,
  `internal-readSampleMetadataFile.R`, `internal-sampleDirs.R`. We may want to
  provide this code as a shared bcbio core package (e.g. [bcbioBase][]) in the
  future.
- Added internal utility to check for valid marker genes (`.validMarkers`).
- Improved [Ensembl][] release version support (`ensemblVersion`).

# bcbioSingleCell 0.0.18 (2017-09-17)

- Renamed `plotClusters` to `plotMarkers`. Added soft deprecation.
- Added [viridis][] color support in t-SNE plots and heatmaps.
- Converted `loadSingleCellRun` and `loadCellRanger` from S4 generics back
  to standard functions.
- Added t-SNE utility functions: `fetchTSNEData`, `fetchTSNEExpressionData`,
  and `plotTSNEExpressionData`. This enable plotting of geometric mean values
  of desired marker genes.
- Updated NEWS to Markdown, with hyperlinks.
- Offloaded generics that would otherwise conflict with bcbioRNASeq to the
  basejump package.
- Improved [roxygen2][] documentation. Moved as much documentation as possible
  to the methods files.
- Updated `cellCycleMarkers` and `cellTypeMarkers` data. Now supports
  Drosophila.
- Sample IDs are now sanitized using `make.names` instead of `camel`. This
  avoids undesirable coercion of some IDs (e.g. `group1_1` into `group11`).
- Added recommended package syntax guidelines.
- lintr checks now allow implicit integers (e.g. `1` instead of `1L`).
- Added [Seurat][] as dependency in `DESCRIPTION` file. The package now attaches
  [Seurat][] automatically.
- Package no longer imports mononcle or suggests scater, scde, or scone. We're
  planning on adding these back in a future update, but build checks on
  [Travis CI][] otherwise take too long.
- Added new `quantileHeatmap` function.
- Improved Markdown header support across functions, where applicable.
- Improved `bcbioSCFiltered` to `seurat` coercion to slot relevant bcbio
  metadata.

# bcbioSingleCell 0.0.17 (2017-09-03)

- Renamed package from `bcbioSinglecell` to `bcbioSingleCell`.
- Added [viridis][] color palette support to quality control plots.
- Added cell-cycle marker genes for human and mouse.
- Slotted `organism` in `bcbioSCDataSet` metadata, in addition to `genomeBuild`.
- Updated `bcbioSCFiltered` to `seurat` coercion method to also run
  `FindVariableGenes` and `ScaleData` by default.
- Added cell-cycle regression into [Seurat][] clustering [RMarkdown][] template.
- Improved [pkgdown][] settings and website appearance.

# bcbioSingleCell 0.0.16 (2017-08-25)

- Support for CRAN release of [Seurat][].
- Improved documentation of package NAMESPACE in `bcbioSinglecell-package.R`
  file.
- Offloaded `download` functionality to [basejump][] package, to avoid
  collisions with [bcbioRNASeq][] package. Function has been renamed to
  `externalFile`.
- Removed function deprecations to simplify the NAMESPACE.
- Improved [Cell Ranger][] sample matching. With these changes the internal
  `.detectPipeline` function is no longer needed.
- Improved sample directory matching in internal `.sampleDirs` function.
- Updated paths to [Cell Ranger][] MatrixMarket files in internal
  `.readSparseCounts` function.
- Updated use of `packageSE` to `prepareSE`, matching the corresponding
  [basejump][] function change.
- Renamed internal use of `filter` to `tidy_filter`, to avoid future
  NAMESPACE collisions with [ensembldb][] package.
- Renamed `bcbioSCSubset` class to `bcbioSCFiltered` class.
- Fixed memory issue in `plotZeroesVsDepth` for datasets with high cell
  counts.
- Improved sample matching for `selectSamples`. Will attempt to migrate this
  to bracket-based subsetting in a future update.
- Restricted sample selection with `selectSamples` to only work on
  `bcbioSCFiltered` class for the time being. We can add bracket-based
  subsetting or S4 method support in `selectSamples` to properly work on
  `bcbioSCDataSet` class in a future update.
- Updated [RMarkdown][] template settings and simplified the setup chunks.
- Initial support for `bcbioSCFiltered` class coercion to [monocle][]
  `CellDataSet` class.
- Suggest [scater][] and [scone][] packages for additional quality control and
  visualization.
- Require [monocle][] >= 2.5.0.
- Added initial support for [viridis][] color palette in quality control plots.
- Suggest [rmarkdown][] and [scde][] packages.
- Initial commit of `subsetPerSample` function.
- Miscellaneous [RMarkdown][] template improvements.

# bcbioSingleCell 0.0.15 (2017-08-11)

- Renamed functions in `lowerCamelCase` from `snake_case`.
- Draft support for [monocle][], [scater][], and [scone][].
- Package now depends on [SummarizedExperiment][].
- Renamed `loadRun` to `loadSingleCellRun` for improved compatibility with
  [bcbioRNASeq][] package. This helps avoid NAMESPACE collisions between
  packages.
- Improved support for custom GTF files.
- Split out [Cell Ranger][] import into a separate utility function named
  `loadCellRanger`.
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

- Integrated [bcbio][] and [Cell Ranger][] workflows into `load_run`.

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
- Changed presentation of mitochondrial abundance as ratio instead of
  percentage.
- Simplified NAMESPACE by offloading dependencies to basejump.

# bcbioSingleCell 0.0.7 (2017-05-13)

- Initial support for loading of 10x Genomics [Cell Ranger][] output.
- Add detection of droplet method based on the metadata file.
- Draft support for import of [inDrop][] i5 index barcode counts.

# bcbioSingleCell 0.0.6 (2017-05-10)

- Modified `load_run` function for improved consistency with [bcbioRNASeq][]
  package.
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

[10x genomics]: https://www.10xgenomics.com/
[acid genomics]: https://acidgenomics.com/
[AcidPlots]: https://AcidPlots.acidgenomics.com/
[basejump]: https://basejump.acidgenomics.com/
[bcbio]: https://bcbio-nextgen.readthedocs.io/
[bcbiobase]: https://bioinformatics.sph.harvard.edu/bcbioBase/
[bcbiornaseq]: https://bioinformatics.sph.harvard.edu/bcbioRNASeq/
[bcbiosinglecell]: https://bioinformatics.sph.harvard.edu/bcbioSingleCell/
[biocinstaller]: https://bioconductor.org/packages/BiocInstaller/
[biocmanager]: https://bioconductor.org/packages/BiocManager/
[bioconductor]: https://bioconductor.org/
[biocparallel]: https://bioconductor.org/packages/BiocParallel/
[cell ranger]: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
[chromium]: https://chromium.acidgenomics.com/
[conda]: https://conda.io/
[covr]: https://github.com/jimhester/covr/
[deseq2]: https://bioconductor.org/packages/DESeq2/
[dplyr]: https://dplyr.tidyverse.org/
[edger]: https://bioconductor.org/packages/edgeR/
[ensembl]: https://www.ensembl.org/
[ensembldb]: https://bioconductor.org/packages/ensembldb/
[ggplot2]: https://ggplot2.tidyverse.org/
[github]: https://github.com/
[google sheets]: https://www.google.com/sheets/
[indrop]: https://1cell-bio.com/
[lintr]: https://github.com/jimhester/lintr/
[macos]: https://www.apple.com/macos/
[monocle]: http://cole-trapnell-lab.github.io/monocle-release/
[pkgdown]: https://pkgdown.r-lib.org/
[pointillism]: https://pointillism.acidgenomics.com/
[r markdown]: http://rmarkdown.rstudio.com/
[r]: https://www.r-project.org/
[rlang]: http://rlang.r-lib.org/
[roxygen2]: https://cran.r-project.org/package=roxygen2
[scater]: https://bioconductor.org/packages/scater/
[scde]: https://bioconductor.org/packages/scde/
[scone]: https://bioconductor.org/packages/scone/
[scran]: https://bioconductor.org/packages/scran/
[seurat]: http://satijalab.org/seurat/
[summarizedexperiment]: https://bioconductor.org/packages/SummarizedExperiment/
[surecell]: https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/surecell-wta-ddseq.html
[testthat]: https://github.com/hadley/testthat/
[tidyeval]: https://dplyr.tidyverse.org/articles/programming.html
[tidyverse]: http://www.tidyverse.org/
[travis ci]: https://travis-ci.org/
[viridis]: https://cran.r-project.org/package=viridis
[zinbwave]: https://bioconductor.org/packages/zinbwave/
[zinger]: https://github.com/statOmics/zingeR/
