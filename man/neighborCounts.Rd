\name{neighborCounts}
\alias{neighborCounts}

\title{Load Hi-C interaction counts}

\description{Collate count combinations for interactions between pairs of bins across multiple Hi-C libraries.}

\usage{
neighborCounts(files, param, width=50000, filter=1L, flank=NULL, 
    exclude=NULL, prior.count=NULL)
}

\arguments{
	\item{files}{a character vector containing paths to the index files generated from each Hi-C library}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{width}{an integer scalar specifying the width of each square in base pairs}
	\item{filter}{an integer scalar specifying the minimum count for each square}
	\item{flank}{an integer scalar, specifying the number of bins to consider as the local neighborhood}
	\item{exclude}{an integer scalar, specifying the number of bins to exclude from the neighborhood}
	\item{prior.count}{a numeric scalar indicating the prior count to use in computing the log-fold increase}
}

\value{
A list object is returned containing \code{interaction}, a \code{DIList} object with the number of read pairs for each bin pair across all libraries; and \code{enrichment}, a numeric vector of enrichment statistics for each bin pair in \code{interaction}.
}

\details{
This function combines the functionality of \code{\link{squareCounts}} and \code{\link{enrichedPairs}}.
The idea is to allow calculation of local enrichment values when there is insufficient memory to load all bin pairs with \code{filter=1L} in \code{\link{squareCounts}}. 
Here, the interaction space around each bin pair is examined as the counts are loaded for that bin pair, avoiding the need to hold the entire interaction space at once.
Only the counts and local enrichment values for those bin pairs with row sums above \code{filter} are reported to save memory.
The returned \code{enrichment} values are equivalent to that computed with \code{\link{enrichedPairs}} with the default settings.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, file=fout))

# Collating to count combinations.
y <- neighborCounts(fout, param, width=50, filter=2, 
    flank=5, prior.count=2)
y$interaction
y$enrichment

# Practically identical to a more memory-intensive call.
ref <- squareCounts(fout, param, width=50, filter=1)
keep <- rowSums(counts(ref)) >= 2
enriched <- enrichedPairs(ref, flank=5, prior.count=2)

stopifnot(identical(ref[keep,], y$interaction))
stopifnot(all(abs(enriched[keep] - y$enrichment) < 1e-6))

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\author{Aaron Lun}

\seealso{
\code{\link{squareCounts}},
\code{\link{enrichedPairs}}
}

\references{
Rao S et al. (2014). A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. \emph{Cell}. 159, 1665-1690.
}

\keyword{counting}
