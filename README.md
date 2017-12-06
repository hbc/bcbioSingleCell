# bcbioSingleCell

[![Build Status](https://travis-ci.org/hbc/bcbioSingleCell.svg?branch=master)](https://travis-ci.org/hbc/bcbioSingleCell)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![codecov](https://codecov.io/gh/hbc/bcbioSingleCell/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioSingleCell)

Import and analyze [bcbio][] single-cell RNA-seq data.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
biocLite("ensembldb")
biocLite(
    "hbc/bcbioSingleCell",
    dependencies = c("Depends", "Imports", "Suggests")
)
```


## Sample metadata examples

### FASTQ files with samples multiplexed by index barcode

This is our current method for handling inDrop samples.

| fileName           | description | index | sequence | sampleName |
| -------------------|-------------|-------|----------|------------|
| 170201_R1.fastq.gz | run1        | 17    | GGAGGTAA | sample1    |
| 170201_R1.fastq.gz | run1        | 18    | CATAACTG | sample2    |
| 170620_R1.fastq.gz | run2        | 12    | GCGTAAGA | sample3    |
| 170620_R1.fastq.gz | run2        | 13    | CTATTAAG | sample4    |
| 170620_R1.fastq.gz | run2        | 14    | AAGGCTAT | sample5    |
| 170620_R1.fastq.gz | run2        | 15    | GAGCCTTA | sample6    |
| 170620_R1.fastq.gz | run2        | 16    | TTATGCGA | sample7    |

### FASTQ files demultiplexed per sample

This is our current method for handling 10X and SureCell samples.

| fileName                   | description |
|----------------------------|-------------|
| sample1_170620_R1.fastq.gz | sample1     |
| sample2_170620_R1.fastq.gz | sample2     |
| sample3_170620_R1.fastq.gz | sample3     |
| sample4_170620_R1.fastq.gz | sample4     |
| sample5_170620_R1.fastq.gz | sample5     |

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

```r
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
[devtools]: https://cran.r-project.org/package=devtools
[Paperpile]: https://paperpile.com
[R]: https://www.r-project.org
