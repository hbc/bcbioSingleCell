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
biocLite("hbc/bcbioSingleCell")
```

### [devtools][] method

```r
install.packages("devtools")
devtools::install_github("hbc/bcbioSingleCell")
```


## Sample metadata examples

### FASTQ files with samples multiplexed by index barcode

This is our current method for handling inDrop samples.

```
fileName	description	index	sequence	sampleName
170201_R1.fastq.gz	run1	17	GGAGGTAA	sample1
170201_R1.fastq.gz	run1	18	CATAACTG	sample2
170620_R1.fastq.gz	run2	12	GCGTAAGA	sample3
170620_R1.fastq.gz	run2	13	CTATTAAG	sample4
170620_R1.fastq.gz	run2	14	AAGGCTAT	sample5
170620_R1.fastq.gz	run2	15	GAGCCTTA	sample6
170620_R1.fastq.gz	run2	16	TTATGCGA	sample7
```

### FASTQ files demultiplexed per sample

This is our current method for handling 10X and SureCell samples.

```
fileName	description	index	sequence
sample1_170620_R1.fastq.gz	sample1	12	GCGTAAGA
sample2_170620_R1.fastq.gz	sample2	13	CTATTAAG
sample3_170620_R1.fastq.gz	sample3	14	AAGGCTAT
sample4_170620_R1.fastq.gz	sample4	15	GAGCCTTA
sample5_170620_R1.fastq.gz	sample5	16	TTATGCGA
```



[bcbio]: https://bcbio-nextgen.readthedocs.io
[bioconductor]: https://bioconductor.org
[devtools]: https://cran.r-project.org/package=devtools
[r]: https://www.r-project.org
