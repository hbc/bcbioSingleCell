globalVariables(".")

packageVersion <- packageVersion("bcbioSingleCell")

# Trailing number is to match cellranger output
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"

requiredAssays <- "counts"
