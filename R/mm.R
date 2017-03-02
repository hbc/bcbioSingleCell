##' Load a gzipped MatrixMart file, setting the row and column names.
##'
##' @param gzfile gzipped MatrixMart file to load
##' @return a sparse matrix of counts
##' @importFrom Matrix readMM
##' @importFrom tools file_path_sans_ext
##' @importFrom readr read_lines
##' @author Rory Kirchner
##' @export
load_sparsecounts = function(gzfile) {
  stem = file_path_sans_ext(gzfile)
  rowfile = paste0(stem, ".rownames")
  colfile = paste0(stem, ".colnames")
  counts = readMM(gzfile)
  rownames(counts) = read_lines(rowfile)
  colnames(counts) = read_lines(colfile)
  return(counts)
}

##' Write a gzipped MatrixMart file, storing the row and column names.
##'
##' @param sparse sparse matrix of counts
##' @param gzfile outfile with gzip extension
##' @return gzfile
##' @importFrom Matrix writeMM
##' @importFrom tools file_path_sans_ext
##' @importFrom R.utils gzip
##' @importFrom readr write_lines
##' @author Rory Kirchner
##' @export
write_sparsecounts = function(sparse, gzfile) {
  if(file.exists(gzfile)) {
    return(gzfile)
  }
  stem = file_path_sans_ext(gzfile)
  rows = rownames(sparse)
  cols = colnames(sparse)
  rowfile = paste0(stem, ".rownames")
  colfile = paste0(stem, ".colnames")
  writeMM(sparse, stem)
  gzip(stem)
  write_lines(rows, rowfile)
  write_lines(cols, colfile)
  return(gzfile)
}
