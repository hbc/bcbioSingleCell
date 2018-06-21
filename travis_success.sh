#!/bin/bash

Rscript -e 'lintr::lint_package()'
Rscript -e 'covr::codecov()'
R CMD BiocCheck .
