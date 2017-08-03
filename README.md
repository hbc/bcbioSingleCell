# bcbioSinglecell

[![Build Status](https://travis-ci.org/hbc/bcbioSinglecell.svg?branch=master)](https://travis-ci.org/hbc/bcbioSinglecell)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![codecov](https://codecov.io/gh/hbc/bcbioSinglecell/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioSinglecell)

Import and analyze [bcbio][] single-cell RNA-seq data.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
biocLite("hbc/bcbioSinglecell")
```

### [devtools][] method

```r
install.packages("devtools")
devtools::install_github("hbc/bcbioSinglecell")
```


## Contribute

- For all changes, fork or create a new branch, then issue a pull request that will be reviewed.
- Do not commit changes directly to `master` branch.
- Support is only provided for the current release version.

### Required checks

```r
lintr::lint_package()
devtools::document()
devtools::check()
BiocCheck::BiocCheck(getwd())
pkgdown::build_site()
```


[bcbio]: https://bcbio-nextgen.readthedocs.io
[bioconductor]: https://bioconductor.org
[devtools]: https://cran.r-project.org/package=devtools
[r]: https://www.r-project.org
