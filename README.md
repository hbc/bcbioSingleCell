# bcbioSingleCell

[![Travis CI](https://travis-ci.org/hbc/bcbioSingleCell.svg?branch=master)](https://travis-ci.org/hbc/bcbioSingleCell)
[![Codecov](https://codecov.io/gh/hbc/bcbioSingleCell/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioSingleCell)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

Import and analyze [bcbio][] single-cell RNA-seq data.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite(
    "hbc/bcbioSingleCell",
    dependencies = c("Depends", "Imports", "Suggests")
)
```


## Load [bcbio][] run

```r
library(bcbioSingleCell)
bcb <- loadSingleCell(
    uploadDir = "bcbio_indrop/final",
    interestingGroups = c("genotype", "treatment"),
    sampleMetadataFile = "sample_metadata.csv",
    organism = "Homo sapiens",
    ensemblVersion = 90L
)
# Back up all data inside bcbioSingleCell object
flatFiles <- flatFiles(bcb)
saveData(bcb, flatFiles)
```

This will return a `bcbioSingleCell` object, which is an extension of the [Bioconductor][] [SingleCellExperiment][SCE] container class.

Parameters:

- `uploadDir`: Path to the [bcbio][] final upload directory.
- `interestingGroups`: Character vector of the column names of interest in the sample metadata, which is stored in the `sampleData()` accessor slot of the `bcbioSingleCell` object. These values should be formatted in camelCase, and can be reassigned in the object after creation (e.g. `interestingGroups(bcb) <- c("batch", "age")`). They are used for data visualization in the quality control utility functions.
- `organism`: Organism name. Use the full latin name (e.g. "Homo sapiens").

Consult `help("loadSingleCell", "bcbioSingleCell")` for additional documentation.


## Sample metadata examples

### FASTQ files with samples multiplexed by index barcode

This is our current method for handling inDrop samples.

| fileName           | description | index | sequence | sampleName |
| -------------------|-------------|-------|----------|------------|
| 170201_R1.fastq.gz | run_1       | 17    | GGAGGTAA | sample_1   |
| 170201_R1.fastq.gz | run_1       | 18    | CATAACTG | sample_2   |
| 170620_R1.fastq.gz | run_2       | 12    | GCGTAAGA | sample_3   |
| 170620_R1.fastq.gz | run_2       | 13    | CTATTAAG | sample_4   |
| 170620_R1.fastq.gz | run_2       | 14    | AAGGCTAT | sample_5   |
| 170620_R1.fastq.gz | run_2       | 15    | GAGCCTTA | sample_6   |
| 170620_R1.fastq.gz | run_2       | 16    | TTATGCGA | sample_7   |

### FASTQ files demultiplexed per sample

This is our current method for handling 10X and SureCell samples.

| fileName              | description |
|-----------------------|-------------|
| S1_170620_R1.fastq.gz | sample_1    |
| S2_170620_R1.fastq.gz | sample_2    |
| S3_170620_R1.fastq.gz | sample_3    |
| S4_170620_R1.fastq.gz | sample_4    |
| S5_170620_R1.fastq.gz | sample_5    |

### Technical replicates

Use `sampleNameAggregate` to assign groupings for technical replicates:

| fileName                  | description   | sampleNameAggregate |
|---------------------------|---------------|---------------------|
| wildtype_L001_R1.fastq.gz | wildtype_L001 | wildtype            |
| wildtype_L002_R1.fastq.gz | wildtype_L002 | wildtype            |
| wildtype_L003_R1.fastq.gz | wildtype_L003 | wildtype            |
| wildtype_L004_R1.fastq.gz | wildtype_L004 | wildtype            |
| mutant_L001_R1.fastq.gz   | mutant_L001   | mutant              |
| mutant_L002_R1.fastq.gz   | mutant_L002   | mutant              |
| mutant_L003_R1.fastq.gz   | mutant_L003   | mutant              |
| mutant_L004_R1.fastq.gz   | mutant_L004   | mutant              |


## Troubleshooting

### Maximal number of DLLs reached

```
Error: package or namespace load failed for 'bcbioSingleCell' in dyn.load(file, DLLpath = DLLpath, ...):
  maximal number of DLLs reached...
```

Depending on your operating system, you may encounter this error about hitting the DLL limit in [R][]. This issue is becoming more common as RNA-seq analysis packages grow increasingly complex. Luckily, we can configure [R][] to increase the DLL limit. Append this line to your `~/.Renviron` file:

```
R_MAX_NUM_DLLS=150
```

For more information on this issue, consult `help("dyn.load")` in the [R][] documentation. The number of loaded DLLs in an [R][] session can be obtained with `getLoadedDLLs()`.


## References

The papers and software cited in our workflows are available as a [shared library](https://paperpile.com/shared/C8EMxl) on [Paperpile][].


[bcbio]: https://bcbio-nextgen.readthedocs.io
[Bioconductor]: https://bioconductor.org
[conda]: https://conda.io
[devtools]: https://cran.r-project.org/package=devtools
[Paperpile]: https://paperpile.com
[R]: https://www.r-project.org
[SCE]: https://doi.org/doi:10.18129/B9.bioc.SingleCellExperiment
