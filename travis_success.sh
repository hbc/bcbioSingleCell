#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    Rscript -e 'lintr::lint_package()'
    Rscript -e 'covr::codecov()'
    R CMD BiocCheck .
fi
