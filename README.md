# bcbioSingleCell

[![Travis CI](https://travis-ci.org/hbc/bcbioSingleCell.svg?branch=master)](https://travis-ci.org/hbc/bcbioSingleCell)
[![Codecov](https://codecov.io/gh/hbc/bcbioSingleCell/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioSingleCell)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[R][] package for [bcbio][] single-cell RNA-seq analysis.

## Installation

This is an [R][] package.

### [Bioconductor][]

We recommend installing the package with [BiocManager][].

```r
if (!require("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install(
    pkgs = c(
        "devtools",
        "remotes",
        "GenomeInfoDbData"
    )
)
BiocManager::install("hbc/bcbioSingleCell")
```

For [R][] < 3.5, [BiocManager][] is not supported. Use `BiocInstaller::biocLite()` instead of `BiocManager::install()`. This requires sourcing the legacy [Bioconductor][] `biocLite.R` script.

```r
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
```

## Load [bcbio][] run

```r
library(bcbioSingleCell)
bcb <- bcbioSingleCell(
    uploadDir = "indrops/final",
    interestingGroups = c("genotype", "treatment"),
    sampleMetadataFile = "sample_metadata.csv",
    organism = "Homo sapiens",
    ensemblRelease = 90L
)
# Back up all data inside bcbioSingleCell object
flat <- flatFiles(bcb)
saveData(bcb, flat, dir = "data")
```

This will return a `bcbioSingleCell` object, which is an extension of the [Bioconductor][] [SingleCellExperiment][SCE] container class. Consult the `bcbioSingleCell()` constructor function documentation for detailed information on the supported parameters:

```r
help(topic = "bcbioSingleCell", package = "bcbioSingleCell")
```

## Sample metadata examples

### FASTQ files with samples multiplexed by index barcode

This is our current recommended method for analyzing an [inDrops][] dataset. The sample index barcodes are multiplexed per FASTQ set. For Illumina sequencing data, the raw binary base call (BCL) data must be converted into FASTQs (split into `R1`-`R4` files) using [bcl2fastq][].

The [inDrops][] library version is automatically detected by [bcbio][], but ensure that the sample index sequences provided match the library version when attempting to create a `bcbioSingleCell` object. A current list of [inDrops v3 index barcodes](https://github.com/seqcloud/seqcloud/blob/master/workflows/bcbio/scrnaseq/harvard_indrop_v3/index_barcodes.csv) is available from [seqcloud][].

Consult the [bcbio][] documentation for more information on how to configure an [inDrops][] run prior to loading into [R][] with the `bcbioSingleCell()` function.

| description | index | sequence | sampleName | aggregate |
|-------------|-------|----------|------------|-----------|
| indrops1    | 17    | GGAGGTAA | sample1    | indrops   |
| indrops1    | 18    | CATAACTG | sample2    | indrops   |
| indrops2    | 12    | GCGTAAGA | sample3    | indrops   |
| indrops2    | 13    | CTATTAAG | sample4    | indrops   |
| indrops2    | 14    | AAGGCTAT | sample5    | indrops   |
| indrops2    | 15    | GAGCCTTA | sample6    | indrops   |
| indrops2    | 16    | TTATGCGA | sample7    | indrops   |

### FASTQ files demultiplexed per sample

This is our current method for handling [10X Genomics Cell Ranger][cellranger] output (using `readCellRanger()`) and [Illumina SureCell][surecell] sample data.

| description | genotype |
|-------------|----------|
| sample1     | wildtype |
| sample2     | knockout |
| sample3     | wildtype |
| sample4     | knockout |

## References

The papers and software cited in our workflows are available as a [shared library](https://paperpile.com/shared/C8EMxl) on [Paperpile][].

[bcbio]: https://bcbio-nextgen.readthedocs.io
[bcl2fastq]: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
[Bioconductor]: https://bioconductor.org
[BiocManager]: https://cran.r-project.org/package=BiocManager
[CellRanger]: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
[conda]: https://conda.io
[devtools]: https://cran.r-project.org/package=devtools
[inDrops]: https://github.com/indrops/indrops
[Paperpile]: https://paperpile.com
[R]: https://www.r-project.org
[SCE]: https://doi.org/doi:10.18129/B9.bioc.SingleCellExperiment
[SureCell]: https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/surecell-wta-ddseq.html
[seqcloud]: http://seq.cloud
