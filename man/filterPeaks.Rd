\name{filterPeaks}
\alias{filterPeaks}

\title{Filter bin pairs for likely peaks}
\description{Identify bin pairs that are likely to represent punctate peaks in the interaction space.}

\usage{
filterPeaks(data, enrichment, min.enrich=log2(1.5), min.count=5, min.diag=2L, ...)
}

\arguments{
\item{data}{an InteractionSet object produced by \code{\link{squareCounts}} or \code{\link{neighborCounts}}}
\item{enrichment}{a numeric vector of enrichment values, produced by \code{\link{enrichedPairs}} or \code{\link{neighborCounts}}}
\item{min.enrich}{a numeric scalar indicating the minimum enrichment score for a peak}
\item{min.count}{a numeric scalar indicating the minimum average count for a peak}
\item{min.diag}{an integer scalar specifying the minimum diagonal in the interaction space with which to consider a peak}
\item{...}{other arguments to be passed to \code{\link{aveLogCPM}} for the average count filter}
}

\value{
A logical vector indicating whether or not each bin pair is to be considered as a peak.
}

\details{
Filtering on the local enrichment scores identifies high-intensity islands in the interaction space.
However, this alone is not sufficient to identify sensible peaks.
Filtering on the absolute average counts prevents the calling of low-abundance bin pairs with high enrichment scores due to empty neighborhoods.
Filtering on the diagonals prevents calling of high-abundance short-range interactions that are usually uninteresting.
If either \code{min.count} or \code{min.diag} are \code{NULL}, no filtering will be performed on the average counts and diagonals, respectively.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{squareCounts}},
\code{\link{enrichedPairs}},
\code{\link{neighborCounts}}
}

\examples{
# Setting up the object.
a <- 10
b <- 20
regions <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))

set.seed(23943)
all.anchor1 <- sample(length(regions), 50, replace=TRUE)
all.anchor2 <- as.integer(runif(50, 1, all.anchor1+1))
data <- InteractionSet(matrix(rnbinom(200, mu=10, size=10), 50, 4), 
    GInteractions(anchor1=all.anchor1, anchor2=all.anchor2, regions=regions, mode="reverse"), 
    colData=DataFrame(totals=runif(4, 1e6, 2e6)), metadata=List(width=1))

# Getting peaks.
enrichment <- enrichedPairs(data)
summary(filterPeaks(data, enrichment, min.enrich=0.5))
summary(filterPeaks(data, enrichment, min.enrich=0.5, min.count=10))
summary(filterPeaks(data, enrichment, min.enrich=0.5, min.diag=NULL))
}
