% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/plotReadsPerCell-methods.R
\name{plotReadsPerCell}
\alias{plotReadsPerCell}
\alias{plotReadsPerCell,bcbioSingleCell-method}
\title{Plot read counts per cell}
\usage{
plotReadsPerCell(object, ...)

\S4method{plotReadsPerCell}{bcbioSingleCell}(
  object,
  interestingGroups = NULL,
  geom = c("histogram", "ecdf", "violin", "ridgeline", "boxplot"),
  cutoffLine = FALSE,
  title = "Reads per cell"
)
}
\arguments{
\item{object}{Object.}

\item{interestingGroups}{\code{character}.
Groups of interest to use for visualization.
Corresponds to factors describing the columns of the object.}

\item{geom}{\code{character(1)}.
Plot type. Uses \code{\link[base:match.arg]{match.arg()}} internally and defaults to
the first argument in the \code{character} vector.}

\item{cutoffLine}{\code{logical(1)}.
Include a line marking the cutoff.}

\item{title}{\code{character(1)}.
Title.}

\item{...}{Additional arguments.}
}
\value{
\code{ggplot}.
}
\description{
Plot the distribution of read counts for all unfiltered cellular barcodes.
}
\note{
Updated 2023-12-04.
}
\examples{
data(bcb)

## bcbioSingleCell ====
plotReadsPerCell(bcb, geom = "histogram")
plotReadsPerCell(bcb, geom = "ecdf")
}
\author{
Michael Steinbaugh, Rory Kirchner
}
