# bcbioSingleCell

[![Install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/r-bcbiosinglecell/README.html)
![Lifecycle: retired](https://img.shields.io/badge/lifecycle-retired-red.svg)

**NOTE: [bcbio-nextgen][bcbio] is no longer under active development.**
Refer to the [notice of discontinuation][] for additional details.

[R][] package for [bcbio][] single-cell RNA-seq analysis.

## Installation

This is an R package.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
install.packages(
    pkgs = "bcbioSingleCell",
    repos = c(
        "https://r.acidgenomics.com",
        BiocManager::repositories()
    ),
    dependencies = TRUE
)
```

### [Conda][] method

Configure [Conda][] to use the [Bioconda][] channels.

```sh
# Don't install recipe into base environment.
conda create --name='r-bcbiosinglecell' 'r-bcbiosinglecell'
conda activate 'r-bcbiosinglecell'
R
```

## Load bcbio single-cell RNA-seq data

```r
library(bcbioSingleCell)
object <- bcbioSingleCell(
    uploadDir = file.path("indrops", "final"),
    interestingGroups = c("genotype", "treatment"),
    sampleMetadataFile = "sample_metadata.csv",
    organism = "Homo sapiens",
    ensemblRelease = 90L
)
```

This will return a `bcbioSingleCell` object, which is an extension of the
[Bioconductor][] [SingleCellExperiment][sce] container class. Consult the
`bcbioSingleCell()` constructor function documentation for detailed information
on the supported parameters:

```r
help(topic = "bcbioSingleCell", package = "bcbioSingleCell")
```

## Sample metadata examples

### FASTQ files with samples multiplexed by index barcode

This is our current recommended method for analyzing an inDrops dataset.
The sample index barcodes are multiplexed per FASTQ set. For Illumina
sequencing data, the raw binary base call (BCL) data must be converted into
FASTQs (split into `R1`-`R4` files) using [bcl2fastq][].

The inDrops library version is automatically detected by bcbio, but ensure that
the sample index sequences provided match the library version when attempting to
create a `bcbioSingleCell` object.

Consult the bcbio documentation for more information on how to configure an
inDrops run prior to loading into R with the `bcbioSingleCell()` function.

| description | index | sequence | sampleName | aggregate | genotype |
| ----------- | ----- | -------- | ---------- | --------- | -------- |
| indrops1    | 1     | CTCTCTAT | sample1_1  | sample1   | wildtype |
| indrops1    | 2     | TATCCTCT | sample2_1  | sample2   | knockout |
| indrops1    | 3     | GTAAGGAG | sample3_1  | sample3   | wildtype |
| indrops1    | 4     | ACTGCATA | sample4_1  | sample4   | knockout |
| indrops2    | 1     | CTCTCTAT | sample1_2  | sample1   | wildtype |
| indrops2    | 2     | TATCCTCT | sample1_2  | sample2   | knockout |
| indrops2    | 3     | GTAAGGAG | sample1_2  | sample3   | wildtype |
| indrops2    | 4     | ACTGCATA | sample1_2  | sample4   | knockout |

Note that bcbio currently outputs the reverse complement index sequence in the
sample directory names (e.g. `sample-ATAGAGAG`). Define the forward index
barcode in the `sequence` column here, not the reverse complement. The reverse
complement will be calculated automatically and added as the `revcomp` column
in the sample metadata.

### FASTQ files demultiplexed per sample

This is our current method for handling 10X Genomics Chromium and Illumina
SureCell cell barcodes.

| description | genotype |
| ----------- | -------- |
| sample1     | wildtype |
| sample2     | knockout |
| sample3     | wildtype |
| sample4     | knockout |

### Invalid object

If you encounter a `validObject` error when attempting to load a
`bcbioSingleCell` object from a previous analysis, run this step to update the
object to the current version of the package:

```r
object <- updateObject(object)
validObject(object)
## [1] TRUE
```

## References

The papers and software cited in our workflows are available as a [shared
library](https://paperpile.com/shared/C8EMxl) on [Paperpile][].

[bcbio]: https://bcbio-nextgen.readthedocs.io/
[bcl2fastq]: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
[bioconda]: https://bioconda.github.io/
[bioconductor]: https://bioconductor.org/
[conda]: https://conda.io/
[notice of discontinuation]: https://github.com/bcbio/bcbio-nextgen/issues/3749
[paperpile]: https://paperpile.com/
[r]: https://www.r-project.org/
[sce]: https://bioconductor.org/packages/SingleCellExperiment/
