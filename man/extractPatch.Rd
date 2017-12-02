\name{extractPatch}
\alias{extractPatch}

\title{Extract a patch of the interaction space}
\description{Extract and count read pairs into bin pairs for a subset of the interaction space.}

\usage{
extractPatch(file, param, first.region, second.region=first.region, 
    width=10000, restrict.regions=FALSE)
}

\arguments{
\item{file}{character string specifying the path to an index file produced by \code{\link{preparePairs}}}
\item{param}{a \code{pairParam} object containing read extraction parameters}
\item{first.region}{a \code{GRanges} object of length 1 specifying the first region}
\item{second.region}{a \code{GRanges} object of length 1 specifying the second region}
\item{width}{an integer scalar specifying the width of each bin in base pairs}
\item{restrict.regions}{A logical scalar indicating whether the output regions should be limited to entries in \code{param$restrict}.}
}

\details{
This function behaves much like \code{\link{squareCounts}}, but only for the \dQuote{path} of the interaction space defined by \code{first.region} and \code{second.region}.
Read pairs are only counted into bin pairs where one end overlaps \code{first.region} and the other end overlaps \code{second.region}.
This allows for rapid extraction of particular regions of interest without having to count across the entire interaction space.

Note that the first anchor region (i.e., bin) in each bin pair is not necessarily the bin that overlaps \code{first.region}.
In each pair, the bins are sorted so that the first bin has a higher genomic coordinate than the second bin.
The \code{flipped} flag in the metadata of the output object indicates whether this order is flipped.
If \code{TRUE}, the first bin in each pair corresponds to \code{second.region}, and vice versa.

If \code{restrict.regions=TRUE}, only bins on the chromosomes in \code{first.region} and \code{second.region} will be reported in the \code{regions} slot of the output object.
This avoids the overhead of constructing many bins when only a small subset of them are used.
By default, \code{restrict.regions=FALSE} to ensure that the anchor IDs of the output object are directly comparable between different calls to \code{extractPatch}.
}

\value{
An InteractionSet object containing the number of read pairs for each bin pair in the specified patch.
}

\seealso{
\code{\link{squareCounts}}
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- tempfile(fileext=".h5")
invisible(preparePairs(hic.file, param, file=fout))

stuff <- extractPatch(fout, param, GRanges("chrA:1-100"))
interactions(stuff)

stuff <- extractPatch(fout, param, GRanges("chrA:1-100"), GRanges("chrB:1-20"))
interactions(stuff)
}
