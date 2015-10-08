\name{Filtering diagonals}
\alias{filterDiag}

\title{Filtering of diagonal bin pairs}
\description{Filtering to remove bin pairs on or near the diagonal of the interaction space.}

\usage{
filterDiag(data, by.dist=0, by.diag=0L, dist, ...)
}

\arguments{
\item{data}{a \code{DIList} object produced by \code{\link{squareCounts}}}
\item{by.dist}{a numeric scalar indicating the base-pair distance threshold below which bins are considered local}
\item{by.diag}{an integer scalar indicating the bin distance threshold below which bins are considered local}
\item{dist}{a optional numeric vector containing pre-computed distances}
\item{...}{other arguments to pass to \code{\link{getDistance}}, if \code{dist} is not specified}
}

\details{
Pairs of the same bin will lie on the diagonal of the interaction space.
Counts for these pairs can be affected by local artifacts (e.g., self-circles, dangling ends) that may not have been completely removed during earlier quality control steps.
These pairs are also less interesting, as they capture highly local structure that may be the result of non-specific compaction.
In many cases, these bin pairs are either removed or, at least, normalized separately within the analysis.

This function provides a convenience wrapper in order to separate diagonal bin pairs from those in the rest of the interaction space.
Users can also consider near-diagonal bin pairs, which are defined as pairs of local bins on the linear genome.
Specifically, bins are treated as local if they separated by less than \code{by.dist} in terms of base pairs, or by less than \code{by.diag} in terms of bins.
These can be separated with the diagonal bin pairs if they are subject to the same issues described above.

Note that if \code{by.dist} is specified, it should be set to a value greater than 1.5 times the average bin size.
Otherwise, the distance between the midpoints of adjacent bins will always be larger than \code{by.dist}, such that no near-diagonal bin pairs are removed.

Users can expedite processing by supplying a pre-computed vector of distances in \code{dist}.
This vector may already be available if it was generated elsewhere in the pipeline.
However, the supplied vector should have the same number of entries as that in \code{data}.
}

\value{
A logical vector indicating whether each bin pair in \code{data} is a non-diagonal (or non-near-diagonal) element.
}

\seealso{
\code{\link{getDistance}}
}

\author{
Aaron Lun
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, file=fout))

# Collating to count combinations.
y <- squareCounts(fout, param, width=50, filter=1)

summary(filterDiag(y))
summary(filterDiag(y, by.dist=100))
summary(filterDiag(y, by.diag=1))
summary(filterDiag(y, dist=getDistance(y)))

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\keyword{filtering}
