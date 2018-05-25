# bcbioSingleCell 0.1.12 (2018-05-25)

## Major changes

- No longer using automatic camel case sanitization for `metrics()` or
  `fetchData` return column names.
- Improved R Markdown clustering and marker templates to optionally support
  UMAP and dark mode in the YAML parameters.
  
## Minor changes

- Using original [Seurat][] mapping names for data: tSNE_1, tSNE_2, PC1, PC2,
  UMAP1, UMAP2.
- Ensure transcript-level counts always have `stripTranscriptVersions()` command
  applied, to remove the Ensembl transcript versions if present.
- No longer using labels (e.g. A, B, C, D) on `ggplot` grid return.
- Note that "phase" has been renamed to "Phase" in the R Markdown clustering
  for cell-cycle regression PCA.



# bcbioSingleCell 0.1.11 (2018-05-23)

## New functions

- UMAP is now supported. This functionality is provided in: `plotUMAP()`,
  `plotMarkerUMAP()`, and `plotFeatureUMAP()`. Corresponding fetch functions,
  `fetchUMAPData()` and `fetchUMAPExpressionData()`, have also been added.
- `plotGene()`: Added `seurat` method support. If advanced customization of
  the plot is needed, use `plotDot()` or `plotViolin()` instead, or refer to the
  [Seurat][] documentation for alternates.

## Major changes

- Dimensional reduction and marker plots no longer use dark mode by default.
  The default color palette support for marker plots has been improved to
  consistently use viridis.
- `diffExp()`: improved internal code to work directly on
  `SingleCellExperiment`, removing the need to pass `design` and `group`
  parameters internally. Also added unit testing against [zinbwave][],
  [zingeR][], and [edgeR][] support. [DESeq2][] is supported but runs slowly.
- Reworked `plotFeature()` and `plotMarker()` family of functions. Improved the
  color palette support when `dark = FALSE`, now using a flipped viridis plasma
  color palette.
- `aggregateReplicates()` function has been reworked to return a
  `SingleCellExperiment` object instead of `bcbioSingleCell`. The v0.2.4
  update of [bcbioRNASeq][] behaves similarly with this generic.

## Minor changes

- Reworked the internal handling of some `seurat` `SingleCellExperiment` method
  support, using `as(x, "SingleCellExperiment")` internally, which uses the
  new `Seurat::Convert()` function.
- Made some previously deprecated functions now defunct: `plotClusters()`,
  `plotTSNEExpressionData()`, `loadSingleCellRun()`, `darkTheme()`,
  `pcCutoff()`, `quantileHeatmap()`, `plotKnownMarkers()`, `readMarkers()`,
  `readMarkersFile()`.
- Made `plotFeatures()`, `plotMarker()`, and `plotMarkers()` functions defunct.
- `plotPCElbow()` now returns a plot grid.
- `sanitizeMarkers()`: improved internal code for supported bcbio stashed
  metadata, including `rowRanges`.



# bcbioSingleCell 0.1.10 (2018-05-19)

## Minor changes

- `plotCellTypesPerCluster()` is using `dark = TRUE` by default again.
- Fixed `cell2sample()` handling for multiplexed Cell Ranger data loaded up with
  `readCellRanger()`. Need to use stashed `cell2sample` factor saved in
  `metadata()`, rather than attempting to calculate on the fly with
  `mapCellsToSamples()`.
- Updated [Travis CI][] build checks to include bioc-release on macOS.



# bcbioSingleCell 0.1.9 (2018-05-18)

## Major changes

- No longer attempting to sanitize the rownames for `seurat` objects in coercion
  method. This helps maintain the gene symbol appearance in plotting functions
  for genes with hyphens in the names.

## Minor changes

- Using `BiocParallel::SerialParam()` internally for zinbwave in `diffExp()`.
- Simplified `cell2sample()` internal code to always use `mapCellsToSamples()`
  instead of attempting to use a stashed `vector` inside `metadata()` for
  `SingleCellExperiment` method.
- Removed internal `.applyFilterCutoffs()`, which is no longer necessary since
  this functionality is supported in the S4 subset method.
- Simplified assert checks inside `fetchGene()` functions.
- `plotCellTypesPerCluster()`: revert back to `dark = TRUE` by default.
- Consolidated `plotMarker` and `plotFeature` functions in the documentation.
- `sanitizeMarkers()`: Improved gene identifier matching.
- `topMarkers()` now defaults to `coding = FALSE` by default, since not all
  datasets will contain biotype information.



# bcbioSingleCell 0.1.8 (2018-05-16)

## Minor changes

- Initial `updateObject()` method support for `bcbioSingleCell` class.
- Relaxed `validObject()` validity check to not require sample-level metadata in
  `colData()` yet.



# bcbioSingleCell 0.1.7 (2018-05-15)

## New functions

- Added support for Uniform Manifold Approximation and Projection (UMAP) with
  the `plotUMAP()` and `fetchUMAPData()` functions. These work similarly to the
  other `plotDimensionalReduction()` and `fetchData()` functions.

## Major changes

- Now adding sample-level metadata into `colData()` slot, for better downstream
  compatibility with other packages that work with `SingleCellExperiment`
  container class. Unique per-sample rows are still saved internally in the
  `sampleData()` slot.
- Now recommending ECDF as the default geom for quality control plots, where
  applicable.
- `filterCells()` now supports `minUMIs = c("knee", "inflection")` for automatic
  filtering based on the cellular barcode ranks. Internally this is handled
  by `DropletUtils::barcodeRanks()`.
  
## Minor changes

- Attempting to re-enable `libgsl-dev` installation for [zinbwave][] on
  [Travis CI][].
- Suggesting [BiocParallel][] for [zinbwave][] call in `diffExp()`.
- Now importing [Seurat][] functions into NAMESPACE.
- Consolidated `fetchData()` functions in the documentation.
- Consolidated `plotDimensionalReduction()` functions in the documentation.
- Updated `aggregateReplicates()` internal code. This function again only
  supports aggregation of `bcbioSingleCell` objects that have been filtered
  using the `filterCells()` function.
- Now using `Seurat::Convert()` internally to coerce `seurat` class object
  to `SingleCellExperiment`, using `as(seurat, "SingleCellExperiment")`. This
  utility function was added to [Seurat][] v2.3.1.
  
# Internal changes

- Tweaked `metrics()` `SingleCellExperiment` method code to always merge
  `colData()` and `sampleData()`.
- Updated `readCellRanger()` internal code to match `bcbioSingleCell()`
  constructor, specifically handling sample-level metadata in `colData()`.



# bcbioSingleCell 0.1.6 (2018-05-09)

## Minor changes

- Updated default QC R Markdown template.
- Added trendline option to QC scatterplot functions.
- Simplified internal handling of interestingGroups in `plotQC()`.
- Using `readYAMLSampleData()` internally instead of defunct `sampleYAMLMetadata()`.
- Added `sampleNames()` method support for seurat.
- Now importing rmarkdown, sessioninfo, tidyverse for R Markdown reports, rather
  than suggesting. Similar update applied to bcbioRNASeq.
- Suggesting scater and scran.



# bcbioSingleCell 0.1.5 (2018-05-04)

## Major changes

- Overhauled inflection and knee point labeling support in `plotUMIsPerCell()`.
  Now uses the `point` argument and always labels per sample. Currently requires
  the `geom = "ecdf"` argument for labeling.
- Updated default quality control template.
- Added `plotBarcodeRanks()`.

## Minor changes

- Added barcode rank support for `seurat` class objects.
- QC plots now have titles by default, matching the conventions used in bcbioRNASeq.
- Fixed y-axis scale for histogram geom in QC plots.
- Prefiltering of very low quality barcodes with no UMIs or genes is now always applied. This helps avoid unwanted downstream errors with zero count barcodes.
- Added boxplot geom support for `plotReadsPerCell().
- `plotQC()` geom argument is now more consistent across the paneled plots.
- Fixed facet wrapping for aggregate samples in the QC plots.
- Added `interestingGroups` support to `plotZerosVsDepth()`, matching the other QC functions.



# bcbioSingleCell 0.1.4 (2018-04-30)

- Updated `sampleData()` S4 methods to match update in bcbioBase. Now supports
  `clean` argument, which returns non-blacklisted factor columns only. See
  `bcbioBase::metadataBlacklist` for the blacklist.
- Improved axis scale appearance on dimensionality reduction plots using
  `scales::pretty_breaks()` internally.
- Added `grid` argument to plots, where applicable.
- Renamed example dataset from `bcb_small` to `indrops_small`.
- Removed unnecessary method support for `interestingGroups()` and `metadata()`.
  These extend from `SummarizedExperiment` correctly now.
- Fixed x-axis label centering for `plotCellCounts()` and `plotReadsPerCell()`.
- Simplified seurat method support for `SingleCellExperiment`-like methods,
  where applicable. This includes `rowData`, `gene2symbol()`, and
  `interestingGroups()`.
- Improved dark mode color support for `plotDot()`, `plotFeatureTSNE()`,
  `plotMarker()`.
- Updated sample metadata example in README.



# bcbioSingleCell 0.1.3 (2018-04-25)

## Minor changes

- Improved summary statistics output during `filterCells()` call.
- Miscellaneous documentation improvements, most notably to `bcbioSingleCell()`
  constructor function.
- `plotViolin()` now uses a color border by default.
- Improved `cell2sample` mapping internally for `readCellRanger()`.
- Improved [Bioconductor][] 3.7 installation instructions.



# bcbioSingleCell 0.1.2 (2018-04-24)

## Major changes

- Now using `bcbioSingleCell()` instead of `loadSingleCell()` as the main
  constructor function to create a `bcbioSingleCell` object. `loadSingleCell()`
  is deprecated and still works, but will warn the user.
- Renamed `loadCellRanger()` to `readCellRanger()` for better name consistency.
- Quality control function color palettes now default to [ggplot2][] colors
  instead of using [viridis][] palettes. This is defined using
  `scale_color_hue()` instead of `scale_color_viridis()` for example. The
  [viridis][] color palette is still used by default for marker expression
  plots.
- Use "`aggregate`" instead of "`sampleNameAggregate`" to define
  aggregate/grouped samples in metadata.

## Minor changes

- Reexporting relevant [ggplot2][] and [viridis][] color palettes.
- Renamed references to "inDrops" from "inDrop", where applicable.
- Consistently using `sym()` in place of `.data` internally for tidy code.
- Updated `seurat` blacklist for `sampleData()` generic.
- `plotCellTypesPerCluster()` and `plotMarkerTSNE()` now use an automatic color
  palette by default, which enables for dynamic color palette support when
  `dark = TRUE`. Internally this is handled with the `theme_midnight()` and
  `theme_paperwhite()` [ggplot2][] themes.
- Updated installation instructions to support [Bioconductor][] 3.7.

## Internal changes

- `metrics()` method support now defaults to `matrix` and works similarly
  for `dgCMatrix` sparse matrices. This is used in place of `calculateMetrics()`
  to generate the per cell quality control metrics.
- Reworked and improved `aggregateReplicates()` internal code.
- Fixed facet wrapping when `aggregate` is defined in metadata for quality
  control plots.
- Improved internal code for `sanitizeMarkers()` to use map the gene annotations
  from `rowRanges()` better.



# bcbioSingleCell 0.1.1 (2018-04-16)

## Major changes

- Added support for calculating `barcodeRanks()` and `barcodeRanksPerSample()`.
- Now exporting `plotMarker()` in addition to `plotMarkers()`.
- Primary counts matrix slotted into `assay()` is named `counts` instead of
  `raw`, for better consistency with `SingleCellExperiment` class. The
  `counts()` generic requires that the primary assay slot is named `counts` to
  work correctly. Nothing else here has changed, just the name.
- `loadCellRanger()` now returns a `SingleCellExperiment` object instead of a
  `bcbioSingleCell` object.
- Added support for `transgeneNames` and `spikeNames` when loading up a dataset.
- Datasets from a poorly annotated genome can now be loaded up using
  `organism = NULL` during the `loadSingleCell()` call.
- Methods now dispatch on `SingleCellExperiment` rather than `bcbioSingleCell`
  where applicable, providing support for `SingleCellExperiment` objects created
  elsewhere.

## Minor changes

- Added support for labeling barcode ranks, such as elbow or inflection point on
  UMI counts per cell plots.
- Added support for plotting metrics using an empirical distribution function
  (ECDF) plot.
- Use `theme_midnight()` and `theme_paperwhite()` internally for dimensionality
  reduction plots.
- TSNE and PCA plots now use an aspect ratio of 1 by default.
- Improved imports from Matrix, S4Vectors, and methods.
- Switched back to using base `stop()`, `warning()`, and `message()`.
- `inflectionPoint()` has been made defunct, in favor of using `barcodeRanks()`.
- Moved internal constructors into S4 methods, where applicable.



# bcbioSingleCell 0.1.0 (2018-04-04)

## Major changes

- `bcbioSingleCell` S4 class now extends `SingleCellExperiment` instead of
  `SummarizedExperiment`. This requires definition of `rowRanges()` inside the
  object instead of `rowData()`. Similar functionality was added to the
  bcbioRNASeq package. Upgrade support will be provided using `updateObject()`in
  a future release.
- Added a differential expression utility function named `diffExp()`, which uses
  [zingeR][]/[edgeR][] internally to calculate gene expression changes across
  cell groups.

## Minor changes

- Added `plotCumulativeUMIsPerCell()` utility. This may be removed in a future
  update in favor of adding this plot into `plotUMIsPerCell()` using an ECDF
  plot.
- Now using `readCellTypeMarkers()` to load marker `data.frames`, instead of
  `readCellTypeMarkersFile()`. This matches the conventions used in the
  [bcbioBase][] package.



# bcbioSingleCell 0.0.32 (2018-03-12)

- Fixes for object subsetting and `loadSingleCell()` organism calls



# bcbioSingleCell 0.0.31 (2018-02-21)

- Improved handling of old [Ensembl][] release version for [Cell Ranger][]
  output.
- Switched [R Markdown][] templates to consistent snake case formatting for
  variables and file names.
- `prepareSingleCellTemplate()` now uses `_setup.R` instead of `setup.R`.
- Updated unit tests to work with new assert checks.
- Export `mapCellsToSamples()` instead of `cell2sample()`. `cell2sample()` now
  simply acts as an accessor function, returning the internally stored
  cell2sample mappings rather than trying to calculate. `mapCellsToSamples()`
  performs the actual mapping from cellular barcodes to sample identifiers.
- Export `metricsPerSample()`
- Now using `BiocParallel::bpmapply()` to loop across the sparse matrix files
  per sample internally in the `.sparseCountsList()` function, which is shared
  between `loadSingleCell()` and `loadCellRanger()`.
- Updated assert checks.
- Updated `prepareSingleCellTemplate()` to explicitly state which files to
  include for each R Markdown template, rather than inheriting from
  `bcbioBase::prepareTemplate()`.
- `selectSamples()` now fails on a sample mismatch, rather than warning.
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
- Added new `fetchGeneData()` function, that wraps the functionality of
  `Seurat::FetchData()` for specific genes.
- Improved internal dimensionality reduction plotting code. The main function
  has been renamed to `.plotDR()` internally.
- Simplified the code for `fetchTSNEExpressionData()` to return a standard
  `data.frame` with the cellular barcodes as rows, instead of the previously
  grouped `tibble` method. Now this function returns aggregate gene marker
  calculations in the `mean`, `median`, and `sum` columns. Since this method has
  way fewer rows than the grouped `tibble`, the [ggplot2][] code for
  `plotMarkerTSNE()` now runs faster.
- Explicitly declare `viridis::` for color palettes, where applicable.
- `colorPoints` argument has been renamed to `expression` for
  `plotMarkerTSNE()`.
- Dynamic gene symbol to ensgene conversion has been removed from the plotting
  functions, for greater simplicity. Now the `genes` argument simply matches
  against the rownames in the counts matrix of the object.
- Decreased the default `minCumPct` argument for `plotPCElbow()` from 0.9 to
  0.8. This is more conservative and will return slightly fewer principal
  components for dimensionality reduction, by default.
- `plotTSNE()`, `plotPCA()`, and the other dimensionality reduction-related
  plotting functions now default to a smaller point size (0.5) and slight alpha
  transparency (0.8), to make super imposed points more obvious for large
  datasets with many cells.
- Added a working example for `subsetPerSample()`.
- Internally switched from `.onLoad()` to `.onAttach()` method for automatically
  loading required dependency packages.
- `plotPCA()` now uses `phase` instead of `Phase` plotting cell cycle regression
  as an interestingGroup (see Seurat clustering template). Previously some of
  the [Seurat][] metadata columns were not consistently sanitized to
  lowerCamelCase (e.g. `Phase`, `res.0.8`, `orig.ident`).
- Suppress package startup messages in [R Markdown][] templates.



# bcbioSingleCell 0.0.29 (2018-01-24)

- Switched to [rlang][] methods for errors, messages, and warnings: `abort()`,
  `inform()`, and `warn()`.
- Updated `filterCells()` function to enable per sample filtering cutoffs. This
  works by passing in a named numeric vector, where the names must match the
  internal `sampleID` metadata column (not `sampleName`).
- Improved internal sanitization of metrics available with the `metrics()`
  accessor. Now all count columns (e.g. `nUMI`) are consistently integers, and
  all character vector columns are consistently coerced to factors.
- Seurat metadata available through the bcbioSingleCell generics are now
  consistently sanitized in lowerCamelCase. This applies to `orig.ident` and the
  `res.*` metadata columns.
- Explicit integers are now consistently used in all of the function parameter
  arguments.



# bcbioSingleCell 0.0.28 (2018-01-22)

- Manually define functions used to read barcodes and matrices. This improves
  the functionality of the internal `.readSparseCounts()` function.
- Improved sample directory matching for [Cell Ranger][] output.
- Fixed sample metadata subsetting based on cell2sample factor levels in the
  internal subset code.
- Added tabbed histogram, violin, and barplots for appropriate quality control
  functions in the [R Markdown][] code.



# bcbioSingleCell 0.0.27 (2018-01-18)

- Migrated core dependency imports from [basejump][] to [bcbioBase][].
- Improved `colnames` and `rownames` handling for internal `.readSparseCounts()`
  function.
- Reworked `loadCellRanger()` function. `refDataDir` parameter has been renamed
  to `refdataDir`.
- Added `organism` and `genomeBuild` options to `loadSingleCell()`, to override
  the metadata set in the bcbio run, if necessary.
- Improved if statement data class checks, where applicable.
- Renamed internal `.sparseCountsTx2Gene()` to `.transcriptToGeneLevelCounts()`.
- Updated `geomean` bind method in `fetchTSNEExpressionData()`.
- Changed `minNovelty` default from 0.8 to 0.75.
- Improved seurat class support for `plotDot()`.
- Improved parameter names for `plotKnownMarkersDetected()`. Now uses
  `tsneColor`, `violinFill`, and `dotColor`. Also added `pointsAsNumbers`
  parameter.
- Added `subtitle` parameter for `plotMarkerTSNE()`.
- Added `tsneColor`, `violinFill`, `dotColor`, and `dark` parameters for
  `plotMarkers()`.
- Improved looping method for `plotTopMarkers()` so that it renders correctly in
  [R Markdown][] calls.
- Added method support for `plotViolin()`.
- Improved internal `cell2sample` handling in subset method code.



# bcbioSingleCell 0.0.26 (2017-12-18)

- Renamed `readMarkersFile()` to `readCellTypeMarkersFile()`.
- Improved internal handling of multiplexed CellRanger samples. These are count
  matrices with cellular barcodes ending in `-2`, for example.
- Updated `aggregateReplicates()` code to work with basejump generic, which uses
  `groupings` instead of `cells` as the grouping parameter.
- Added method support for `detectOrganism()`.
- Improved filtering parameter output in `filterCells()`.
- Updated `gene2symbol()` method support for bcbioSingleCell and seurat objects.
- Fixed working example in `knownMarkersDetected()`.
- Reworked internal code for `plotCellTypesPerCluster()`.
- Improved internal checks for facet wrapping and color palette parameter
  arguments, where applicable.
- Offloaded base `plotQuantileHeatmap()` functionality into basejump, for use in
  [bcbioRNASeq][] package.
- Added unit test data for `loadSingleCell().
- Added additional unit tests to improve code coverage.



# bcbioSingleCell 0.0.25 (2017-12-11)

- Prepared a minimal working example dataset.
- Added some initial unit tests.
- Renamed `plotKnownMarkers()` to `plotKnownMarkersDetected()`.
- Renamed `readMarkers()` to `readMarkersFile()`.
- Added `gene2symbol()` method support.
- Moved `plotDot()` generic to basejump package.
- Added internal `.checkFormat()` function, which will check for `ensgene` or
  `symbol` input.
- Added internal `.convertGenesToSymbols()` utility function, for mapping
  Ensembl gene identifiers to gene symbols.
- Use explicit calls to Seurat functions interally, for clarity.
- Simplified bcbioSingleCell object return in `loadSingleCell()` function.
- Added seurat class method support for `annotable()` function.
- Added `bcbio<-` assignment method support for seurat class objects.
- `calculateMetrics()` function now uses `annotable = TRUE` as default, instead
  of using `missing()` method.
- Added method support for `cell2sample()` for seurat class objects.
- `cellTypesPerCluster()` now uses `min` and `max` default arguments that don't
  remove any rows.
- Added method support for seurat class objects to `counts()` function. This
  defaults to returning the raw counts (`normalized = FALSE`), but can also
  return log-normalized (`normalized = TRUE`) and scaled
  (`normalized = "scaled"`) counts.
- Added Ensembl gene identifier mapping support to `fetchTSNEExpressionData()`.
- `filterCells()` now simply works in a destructive manner. We've removed the
  `drop` parameter. The messages displayed to the user during this function call
  have been improved, and now include more statistics on the step where the
  majority of cells are filtered.
- Initial commit of `gene2symbol()` method support for bcbioSingleCell and
  seurat class objects.
- Improved `interestingGroups<-` assignment method support for seurat class
  objects.
- Color palettes now default to viridis instead of inferno palette, where
  applicable.
- Added Ensembl gene identifier support for `plotDot()` function.
- `plotFeatureTSNE()` now uses a plural `features` parameter instead of
  `feature`, which is consistent with the syntax used in the other functions.
- Added Ensembl gene identifier support to `plotMarkerTSNE()`. The `format`
  argument still defaults to "symbol", for consistency with previous behavior.
  However, in the future we recommend that users pass in stable Ensembl gene
  identifiers here if possible, for better reproducibility.
- `plotMarkers()` now supports Ensembl gene identifiers.
- `plotPCElbow()` now silently returns the sequence of principal components
  (PCs) that we recommend to use for dimensionality reduction.
- We're now using "glm" instead of "gam" for `geom_smooth()` plotting, where
  applicable. See the `plotQC()` function code.
- Improved the defaults for `plotQuantileHeatmap()` to enable faster plotting.
  Now the dendrogram calculations are skipped by default, which take a long time
  for large datasets.
- Removed the draft `plotStressGenes()` function for now. Will try to add this
  in a future update.
- Fixed Markdown header handling for `plotTopMarkers()` if `headerLevel = NULL`.
- Simplified `selectSamples()` code to rely upon output of our bracket-based
  subsetting method. See `subset.R` file for more details.
- Improved metadata update in bracket-based subsetting method.
- `subsetPerSample()` function now defaults to saving in the working directory.
- Updated `topBarcodes()` function to rank by `nUMI` instead of `nCount` column,
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
- Renamed `pcCutoff()` to `plotPCElbow()`. The function now returns a PC
  sequence from 1 to the cutoff (e.g. 1:10) instead of just the final PC cutoff
  value. The R Markdown clustering template has been updated to reflect this
  change.
- Renamed `quantileHeatmap()` to `plotQuantileHeatmap()`, for consistency with
  other plotting functions.
- Initial support for customized bracket based subsetting, which now acts upon
  the raw cellular barcode counts stashed in the `object@bcbio` slot.
- Moved `darkTheme()` to [basejump][] package and reworked as `midnightTheme()`,
  with improved colors and axis appearance.
- Added `pointsAsNumbers` parameter to `plotTSNE()` and `plotPCA()` functions,
  to match the functionality in `plotMarkerTSNE()`.
- Overhauled `loadCellRanger()` to support multiplexed [Cell Ranger][] matrix
  output. [Cell Ranger][] adds a numeric suffix to the end of multiplexed
  barcodes (e.g. `AAACCTGGTTTACTCT-1` denotes cellular barcode
  `AAACCTGGTTTACTCT` is assigned to sample `1`).
- Improved `cell2sample` mapping in `aggregateReplicates()` function, which uses
  the `sampleNameAggregate` column in sample metadata to define the aggregate
  sample pairings. The `summarize()` step at line 101 is slow for datasets with
  many samples and should be changed in the future to speed things up.
- Improved internal `cell2sample()` code to handle `NULL` stashed mappings
  better.
- Updated TNSE plotting functions to use `midnightTheme()` instead of
  `darkTheme()`.
- Added user-defined point and label sizes for `plotMarkerTSNE()`.
- Fixed typo in `plotMitoRatio()` where `maxGenes` cutoff was plotted instead of
  `maxMitoRatio`.
- Added `legend` parameter argument to `plotQC()` function. Also improved
  handling of `NULL` return for `plotReadsPerCell()`, which can happen with
  [Cell Ranger][] output.
- Updated facet wrapping in `plotZeroesVsDepth()` to match the behavior in the
  other plotting functions.
- Initial methods support for custom bracket-based subsetting.



# bcbioSingleCell 0.0.23 (2017-11-22)

- Improved facet wrapping of aggregated samples (`sampleNameAggregate` present
  in sample metadata), but removing code support for wrapping by multiplexed
  FASTQ description.
- Simplified handling of `bcbioSingleCell` objects with `filterCells()` applied.
  This information is stored in the `metadata()` slot as 3 variables: (1)
  `filterParams`, numeric vector of the parameters used to define the cell
  filtering cutoffs; (2) `filterCells`, character vector of the cellular barcode
  IDs that passed filtering; (3) `filterGenes`, character vector of the
  [Ensembl][] gene identifiers that have passed filtering, as determined by the
  `minCellsPerGene` parameter.
- For `filterCells()` return, we're now defaulting to a destructive operation,
  where the columns (cells) and rows (genes) of the object are adjusted to match
  the cells and genes that have passed filtering. Currently this can be adjusted
  with the `drop` argument for testing, but should generally be left as
  `drop = TRUE`.
- We're now slotting a `cell2sample` named factor in the `metadata()` slot,
  which makes downstream quality control operations faster. This is generated on
  the fly for previously saved objects that don't have a stashed `cell2sample`.
- Initial commit of `plotQC()` utility function, which plots multiple quality
  control functions for easy visualization. This defaults to output as a cowplot
  grid (`return = "grid"`), but can alternatively be set to return
  [R Markdown][] code (`return = "markdown"`).
- `aggregateReplicates()` operation has been improved to properly slot raw
  cellular barcodes in `object@bcbio$cellularBarcodes`. The `filterCells` vector
  is adjusted, and `sampleMetadata` factors should be properly releveled.
- The `counts()` accessor simply returns the sparse matrix contained in the
  `assay()` slot. The `filterCells` argument has been removed.
- Messages have been added to the `filterCells()` function, to help the user
  determine at which step the majority of cells are being filtered. We're
  keeping a non-destructive option using `drop = FALSE` for the time being, but
  this will likely be removed for improved simplicity in a future update.
- Updated the internal code for `metrics()` to use a simpler join operation on
  the `colData`, `cell2sample` and `sampleMetadata`.
- Updated facet wrap code in quality control plots to not facet multiplexed
  FASTQ descriptions and simply check for `sampleNameAggregate`.
- Improved appearance of `plotReadsPerCell()` labels and legends. Additionally,
  `plotReadsPerCell()` more efficiently handles the stashed values in the
  `nCount` column of `colData`, for faster plotting that having to rely on
  manipulation of the raw `cellularBarcodes` list stashed in
  `object@bcbio$cellularBarcodes`.
- `sampleMetadata()` return is now consistently sanitized for `bcbioSingleCell`
  and `seurat` objects.
- Minor tweaks to quality control template and setup.R files.
- Added `plotFeatureTSNE()` utility function. This improves on
  `Seurat::FeaturePlot()` and enables the user to overlay the cluster
  identifiers on top of the t-SNE plot. `plotFeatures()` is now deprecated in
  favor of this function.
- Improved internal handling of Seurat data in the `.fetchDimDataSeurat()`
  function. This now keeps the cell ID as the rowname.
- Allow the user to define the color palette (`color`), as well as `pointSize`
  and `labelSize` for `plotPCA()` and `plotTSNE()`.
- Factors are now correctly releveled in `cell2sample()` return.
- Improved internal code for `fetchTSNEExpressionData()`.
- Bug fix for `metrics()` accessor not including the cell ID as rownames.
- Clustering template fixes. Now uses `plotFeatureTSNE()` to assess quality
  control metrics on t-SNE.



# bcbioSingleCell 0.0.22 (2017-11-17)

- Now internally stashing a `cell2sample` data.frame, which helps speed up
  operations on cellular barcode metrics calculations for quality control plots.
- Improved support for optional `annotable`, `ensemblVersion`, `gtfFile`, and
  `sampleMetadataFile` arguments in `loadSingleCell()` function.
- Simplified some of the messages shown during sample loading, in an attempt to
  make them clearer and more informative. Also improved messages shown to the
  user during a `filterCells()` function call.
- The `metrics()` function will now look for a stashed `cell2sample`
  `data.frame`, which speeds up operations for quality control plots.
- Improved handling of sample metadata columns as factors. In particular, levels
  should be correctly updated using `droplevels` in a `selectSamples()` call.
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
- Initial support for `plotCellTypesPerCluster()`.
- Initial support for `plotMitoVsCoding()`. I broke this code out from
  `plotMitoRatio()`. We could opt to keep this in `plotMitoRatio` with a
  `geom = "scatterplot"` argument.
- Initial commit of `.applyFilterCutoffs()` internal function, used to subset
  the object to contain only cells and genes that have passed quailty control
  filtering.
- Improved [Seurat][] `FindAllMarkers()` sanitization.
- Made the quality control plots more modular. Now they support multiple geoms,
  including `boxplot`, `histogram`, `ridgeline`, and `violin` (default). Median
  labels are applied with the internal `.medianLabels()` function.
- Updated error message for [Cell Ranger][] directory structure.
- Cellular barcode columns are no longer split with an underscore. Instead, they
  are kept as a single ACGT string. We're now generating `cellID` to `sampleID`
  matching with a different method. In the future, we'll stash a cell2sample
  data.frame inside the object, that makes this operation faster than the
  current `mclapply()` code.
- Updated annotable support in `loadSingleCell()` and `loadCellRanger()`.
- Added support for handling both gene- and transcript-level counts. The updated
  release of the bcbio single-cell pipeline now outputs at gene level.
- Initial quality control plot support for seurat objects.
- Added support for return of only filtered cells and genes with the
  `counts(filterCells = TRUE)` function.
- Added assignment support for `interestingGroups<-`.
- Initial support for sample metadata generation from seurat object.
- Improved internal code for `selectSamples()`.
- Updated `topMarkers()` to match Seurat v2.1 update.
- Improved `bcbioSingleCell` to `seurat` coercion method with `setAs()`.



# bcbioSingleCell 0.0.20 (2017-10-24)

- Upgraded to [basejump][] 0.1.0 and [Seurat][] 2.1.0 dependencies.
- Improved documentation of NAMESPACE imports per function.
- Switched to base grep functions where applicable (`grepl()`, `gsub()`).
- Use GTF in package documentation rather than GFF. Applies primarily to the
  `loadSingleCell()` import function.
- Restrict class support in S4 methods to `bcbioSingleCell`. Legacy
  `bcbioSCDataSet` class can be upgraded to `bcbioSingleCell` class using
  `as(bcb, "bcbioSingleCell")` coercion.
- Use filtered cell output for `metrics()` and quality control functions by
  default.
- Updated the quality control R Markdown to include `filterCells = FALSE` where
  applicable.
- Draft support for aggregated technical replicates in quality control functions
  using `sampleNameAggregate` column in sample metadata. This doesn't change the
  actual counts values. It only applies to visualization in the quality control
  plots currently.
- Miscellaneous [R Markdown][] template updates. Primarily improvements to the
  setup chunk object loading workflow.
- Removed lintr checks from testthat. This is breaking `devtools::test()`.



# bcbioSingleCell 0.0.19 (2017-10-12)

- Renamed main object class from `bcbioSCDataSet` to `bcbioSingleCell`.
- Cell filtering with `filterCells()` will now slot a named logical vector into
  `metadata(object)[["filteredCells"]]`, which will be used to dynamically
  subset the slotted internal `SummarizedExperiment` data. Now that we're using
  this approach, we can return a modified `bcbioSingleCell` object rather than
  defining a separate `bcbioSCFiltered` class.
- Renamed `loadSingleCellRun()` to `loadSingleCell()`, to match [bcbioRNASeq][]
  package.
- Now allowing implicit integers in our function code.
- Added support for plotting technical replicates. This is handled by
`sampleNameAggregate` in the sample metadata.
- Now using ridgeline plots in place of histograms where applicable.
- Travis CI checks take too long when loading SummarizedExperiment. Hopefully
  this will be fixed in the 3.6 release later this month.
- New internal dark theme (`darkTheme()`), based on the Seurat theme.
- Initial commit of `plotDot()` function, based on `Seurat::DotPlot()`.
- Added new t-SNE plots that allow for consistent cluster labeling.
- Providing legacy support for `bcbioSCDataSet` and `bcbioSCFiltered`, which
  will be deprecated in a future release.
- Offloaded some internal code to basejump, for improved consistency with
  [bcbioRNASeq][] package: `internal-projectDir.R`,
  `internal-readSampleMetadataFile.R`, `internal-sampleDirs.R`. We may want to
  provide this code as a shared bcbio core package (e.g. [bcbioBase][]) in the
  future.
- Added internal utility to check for valid marker genes (`.validMarkers()`).
- Improved [Ensembl][] release version support (`ensemblVersion`).



# bcbioSingleCell 0.0.18 (2017-09-17)

- Renamed `plotClusters()` to `plotMarkers()`. Added soft deprecation.
- Added [viridis][] color support in t-SNE plots and heatmaps.
- Converted `loadSingleCellRun()` and `loadCellRanger()` from S4 generics back
  to standard functions.
- Added t-SNE utility functions: `fetchTSNEData()`, `fetchTSNEExpressionData()`,
  and `plotTSNEExpressionData()`. This enable plotting of geometric mean values
  of desired marker genes.
- Updated NEWS to Markdown, with hyperlinks.
- Offloaded generics that would otherwise conflict with bcbioRNASeq to the
  basejump package.
- Improved [roxygen][] documentation. Moved as much documentation as possible to
  the methods files.
- Updated `cellCycleMarkers` and `cellTypeMarkers` data. Now supports
  Drosophila.
- Sample IDs are now sanitized using `make.names()` instead of `camel()`. This
  avoids undesirable coercion of some IDs (e.g. `group1_1` into `group11`).
- Added recommended package syntax guidelines.
- lintr checks now allow implicit integers (e.g. `1` instead of `1L`).
- Added [Seurat][] as dependency in `DESCRIPTION` file. The package now attaches
  [Seurat][] automatically.
- Package no longer imports mononcle or suggests scater, scde, or scone. We're
  planning on adding these back in a future update, but build checks on
  [Travis CI][] otherwise take too long.
- Added new `quantileHeatmap()` function.
- Improved Markdown header support across functions, where applicable.
- Improved `bcbioSCFiltered` to `seurat` coercion to slot relevant bcbio
  metadata.



# bcbioSingleCell 0.0.17 (2017-09-03)

- Renamed package from `bcbioSinglecell` to `bcbioSingleCell`.
- Added [viridis][] color palette support to quality control plots.
- Added cell-cycle marker genes for human and mouse.
- Slotted `organism` in `bcbioSCDataSet` metadata, in addition to `genomeBuild`.
- Updated `bcbioSCFiltered` to `seurat` coercion method to also run
  `FindVariableGenes()` and `ScaleData()` by default.
- Added cell-cycle regression into [Seurat][] clustering [RMarkdown][] template.
- Improved [pkgdown][] settings and website appearance.



# bcbioSingleCell 0.0.16 (2017-08-25)

- Support for CRAN release of [Seurat][].
- Improved documentation of package NAMESPACE in `bcbioSinglecell-package.R`
  file.
- Offloaded `download()` functionality to [basejump][] package, to avoid
  collisions with [bcbioRNASeq][] package. Function has been renamed to
  `externalFile()`.
- Removed function deprecations to simplify the NAMESPACE.
- Improved [Cell Ranger][] sample matching. With these changes the internal
  `.detectPipeline()` function is no longer needed.
- Improved sample directory matching in internal `.sampleDirs()` function.
- Updated paths to [Cell Ranger][] MatrixMarket files in internal
  `.readSparseCounts()` function.
- Updated use of `packageSE()` to `prepareSE()`, matching the corresponding
  [basejump][] function change.
- Renamed internal use of `filter()` to `tidy_filter()`, to avoid future
  NAMESPACE collisions with [ensembldb][] package.
- Renamed `bcbioSCSubset` class to `bcbioSCFiltered` class.
- Fixed memory issue in `plotZeroesVsDepth()` for datasets with high cell
  counts.
- Improved sample matching for `selectSamples()`. Will attempt to migrate this
  to bracket-based subsetting in a future update.
- Restricted sample selection with `selectSamples()` to only work on
  `bcbioSCFiltered` class for the time being. We can add bracket-based
  subsetting or S4 method support in `selectSamples()` to properly work on
  `bcbioSCDataSet` class in a future update.
- Updated [RMarkdown][] template settings and simplified the setup chunks.
- Initial support for `bcbioSCFiltered` class coercion to [monocle][]
  `CellDataSet` class.
- Suggest [scater][] and [scone][] packages for additional quality control and
  visualization.
- Require [monocle][] >= 2.5.0.
- Added initial support for [viridis][] color palette in quality control plots.
- Suggest [rmarkdown][] and [scde][] packages.
- Initial commit of `subsetPerSample()` function.
- Miscellaneous [RMarkdown][] template improvements.



# bcbioSingleCell 0.0.15 (2017-08-11)

- Renamed functions in `lowerCamelCase` from `snake_case`.
- Draft support for [monocle][], [scater][], and [scone][].
- Package now depends on [SummarizedExperiment][].
- Renamed `loadRun()` to `loadSingleCellRun()` for improved compatibility with
  [bcbioRNASeq][] package. This helps avoid NAMESPACE collisions between
  packages.
- Improved support for custom GTF files.
- Split out [Cell Ranger][] import into a separate utility function named
  `loadCellRanger()`.
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
- Changed presentation of mitochondrial abundance as ratio instead of
  percentage.
- Simplified NAMESPACE by offloading dependencies to basejump.



# bcbioSingleCell 0.0.7 (2017-05-13)

- Initial support for loading of 10x Genomics [Cell Ranger][] output.
- Add detection of droplet method based on the metadata file.
- Draft support for import of [inDrop][] i5 index barcode counts.



# bcbioSingleCell 0.0.6 (2017-05-10)

- Modified `load_run()` function for improved consistency with [bcbioRNASeq][]
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



[bcbioBase]: http://bioinformatics.sph.harvard.edu/bcbioBase
[bcbioRNASeq]: http://bioinformatics.sph.harvard.edu/bcbioRNASeq
[Bioconductor]: https://bioconductor.org
[BiocParallel]: https://doi.org/doi:10.18129/B9.bioc.BiocParallel
[Cell Ranger]: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
[inDrop]: https://1cell-bio.com
[R Markdown]: http://rmarkdown.rstudio.com
[tidyverse]: http://www.tidyverse.org
[basejump]: http://steinbaugh.com/basejump
[covr]: https://github.com/jimhester/covr
[Ensembl]: https://www.ensembl.org
[ensembldb]: http://bioconductor.org/packages/release/bioc/html/ensembldb.html
[Google Sheets]: https://www.google.com/sheets
[lintr]: https://github.com/jimhester/lintr
[monocle]: http://cole-trapnell-lab.github.io/monocle-release/
[pkgdown]: https://github.com/hadley/pkgdown
[rlang]: http://rlang.r-lib.org
[roxygen]: https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html
[scater]: http://bioconductor.org/packages/release/bioc/html/scater.html
[scde]: http://bioconductor.org/packages/release/bioc/html/scde.html
[scone]: https://bioconductor.org/packages/release/bioc/html/scone.html
[Seurat]: http://satijalab.org/seurat
[SummarizedExperiment]: https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html
[SureCell]: https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/surecell-wta-ddseq.html
[testthat]: https://github.com/hadley/testthat
[Travis CI]: https://travis-ci.org
[viridis]: https://cran.r-project.org/web/packages/viridis/index.html
[zinbwave]: https://doi.org/doi:10.18129/B9.bioc.zinbwave
[zingeR]: https://github.com/statOmics/zingeR
