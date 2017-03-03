\name{neighborCounts}
\alias{neighborCounts}

\title{Load Hi-C interaction counts}

\description{Collate count combinations for interactions between pairs of bins across multiple Hi-C libraries.}

\usage{
neighborCounts(files, param, width=50000, filter=1L, flank=NULL, exclude=NULL)
}

\arguments{
	\item{files}{a character vector containing paths to the index files generated from each Hi-C library}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{width}{an integer scalar specifying the width of each square in base pairs}
	\item{filter}{an integer scalar specifying the minimum count for each square}
	\item{flank}{an integer scalar, specifying the number of bins to consider as the local neighborhood}
	\item{exclude}{an integer scalar, specifying the number of bins to exclude from the neighborhood}
}

\value{
An InteractionSet object is returned with the number of read pairs for each bin pair across all libraries.
For each bin pair, the number of read pairs in each neighborhood region is also returned in separate \code{assay} fields.
\code{mcols} contains the size of each neighborhood in terms of the number of bin pairs.
}

\details{
This function combines the functionality of \code{\link{squareCounts}} and \code{\link{enrichedPairs}}.
The idea is to allow counting of neighborhoods when there is insufficient memory to load all bin pairs with \code{filter=1L} in \code{\link{squareCounts}}. 
Here, the interaction space around each bin pair is examined as the counts are loaded for that bin pair, avoiding the need to hold the entire interaction space at once.
Only the counts and local enrichment values for those bin pairs with row sums above \code{filter} are reported to save memory.
The returned \code{assay} matrices are equivalent to that computed with \code{\link{enrichedPairs}} with the default settings.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, file=fout))

# Collating to count combinations.
y <- neighborCounts(fout, param, width=50, filter=2, flank=5)
y

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
