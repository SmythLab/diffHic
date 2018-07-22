\name{readMTX2IntSet}
\alias{readMTX2IntSet}

\title{Read InteractionSet from a BED File and a Set of Matrix Market Files}
\description{Read a set of Matrix Market Exchange Format files from disk and create an \code{InteractionSet}.}

\usage{
readMTX2IntSet(mtx, bed, verbose=TRUE)
}

\arguments{
\item{mtx}{character vector containing the names of the Matrix Market Exchange Format (MTX) files to be read.}
\item{bed}{character string giving the name of the BED file to be read.}
\item{verbose}{logical: if \code{TRUE} progress information will be written to standard output as each file is read.}
}

\value{
An \code{InteractionSet} with regions as read from the BED file.
A row will be included for every unique interaction found in any one of the MTX files.
}

\details{
Each MTX file contains the number of reads supporting each Hi-C interaction for one biological sample.
The MTX files are assumed be in Matrix Market Exchange Format (\url{https://math.nist.gov/MatrixMarket/formats.html}), for example as written by the \code{writeMM} function in the \code{Matrix} package.
Files contain three columns, indicating row number, column number and read count respectively.
Counts are assumed to be symmetric and zeros counts are omitted.

The MTX files may optionally contain an Matrix Market header giving the total number of rows and columns for the contact matrix, or the header can be omitted.

Rows and columns correspond to the regions defined in the BED file.
The BED file is assumed to be in standard BED format with no header lines.
Only the first three columns are used, corresponding to chromosome names, start position (in 0-start coordinates) and end position.
Regions shown in diffHic will show start positions in 1-start closed coordinates, being one-plus the start position in the BED file.

The MTX files and the BED file can be either tab-delimited or space-limited.

This function can, in principle, read and merge any number of MTX files.
However it may be more memory efficient to read individual MTX files and subset to interactions of interest before merging into a combined \code{InteractionSet}.
}

\author{Gordon Smyth}

\keyword{read}
