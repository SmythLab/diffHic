\name{readMTX2IntSet}
\alias{readMTX2IntSet}

\title{Create an InteractionSet from a BED file and Matrix Market files}
\description{Read a set of Matrix Market Exchange Format files from disk and create an \linkS4class{InteractionSet} object.}

\usage{
readMTX2IntSet(mtx, bed, as.integer=TRUE)
}

\arguments{
\item{mtx}{A character vector containg paths to Matrix Market Exchange Format (MTX) files.
Each file contains interaction read counts for one sample.}
\item{bed}{String containing the path to the BED file specifying the genomic regions.}
\item{as.integer}{Logical indicating whether the data should be read as integers.
Otherwise values are stored as double-precision numbers.}
}

\value{
An \linkS4class{InteractionSet} object containing interactions between regions from the BED file (in reverse-strict mode, see \code{\link{GInteractions}}).
Each row corresponds to a unique interaction found in any of the MTX files, and contains the read counts across all files.
}

\details{
Each MTX file is assumed to contain read counts for symmetric matrix representing the two-dimensional interaction space.
Each row and column is assumed to correspond to contiguous bins of the genome, with coordinates specified by \code{bed} in standard BED format.
This function will aggregate counts from all files to create an \linkS4class{InteractionSet} object that mimics the output of \code{\link{squareCounts}}.

The \code{width} value in the \code{metadata} of the output InteractionSet is set to the median width of the regions.
The \code{totals} field in the output \code{colData} is also set to be equal to the sum of the counts in each MTX file.
Note that these settings only make sense if the ContactMatrix objects cover binned regions.

This function can, in principle, read and merge any number of MTX files.
However, for large data sets, consider reading each MTX file separately, subsetting it to interactions of interest and then creating the InteractionSet object.
For example, subsetted contact matrices can be used to create an InteractionSet via \code{\link{mergeCMs}}.
}

\author{Gordon Smyth, with modifications by Aaron Lun}

\seealso{
\code{\link{mergeCMs}}
}

\examples{
library(Matrix)
tmp.loc <- tempfile()
dir.create(tmp.loc)

# Mocking up some MTX and BED files.
set.seed(110000)
A <- rsparsematrix(1000, 1000, density=0.1, symmetric=TRUE, 
    rand.x=function(n) round(runif(n, 1, 100)))
A.name <- file.path(tmp.loc, "A.mtx")
writeMM(file=A.name, A)

B <- rsparsematrix(1000, 1000, density=0.1, symmetric=TRUE, 
    rand.x=function(n) round(runif(n, 1, 100)))
B.name <- file.path(tmp.loc, "B.mtx")
writeMM(file=B.name, B)

GR <- GRanges(sample(c("chrA", "chrB", "chrC"), 1000, replace=TRUE),
    IRanges(start=round(runif(1000, 1, 10000)),
        width=round(runif(1000, 50, 500))))
GR <- sort(GR)
bed.name <- file.path(tmp.loc, "regions.bed")
rtracklayer::export.bed(GR, con=bed.name)

# Reading everything in.
iset <- readMTX2IntSet(c(A.name, B.name), bed.name)
iset
}

\keyword{read}
