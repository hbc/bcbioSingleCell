globalVariables(".")

# Trailing number is to match cellranger output
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"

packageVersion <- packageVersion("bcbioSingleCell")
