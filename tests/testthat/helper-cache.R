lst <- AcidDevTools::cacheTestFiles(
    pkg = .pkgName,
    files = "bcbioSingleCell_0.1.0.rds"
)
cacheDir <- lst[["cacheDir"]]
rm(lst)
