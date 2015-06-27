\name{Filtering methods}
\alias{filterDirect}
\alias{filterTrended}

\title{Filtering strategies for bin pairs}
\description{Implementations of the direct and trended filtering strategies for bin pair abundances.}

\usage{
filterDirect(data, prior.count=2, reference=NULL)
filterTrended(data, span=0.25, prior.count=2, reference=NULL)
}

\arguments{
\item{data}{a \code{DIList} object produced by \code{\link{squareCounts}}}
\item{span}{a numeric scalar specifying the bandwidth for loess curve fitting}
\item{prior.count}{a numeric scalar indicating the prior count to use for calculating the average abundance}
\item{reference}{another \code{DIList} object, usually containing data for larger bin pairs}
}

\details{
The \code{filterDirect} function implements the direct filtering strategy.
The rate of non-specific ligation is estimated as the median of average abundances from inter-chromosomal bin pairs.
This rate or some multiple thereof can be used as a minimum threshold for filtering, to keep only high-abundance bin pairs.
When calculating the median, some finesse is required to consider empty parts of the interaction space, i.e., areas that are not represented by bin pairs.

The \code{filterTrended} function implements the trended filtering strategy.
The rate of non-specific compaction is estimated by fitting a trend to the average abundances against the log-distance for all intra-chromosomal bin pairs.
This rate can then be used as a minimum threshold for filtering.
For inter-chromosomal bin pairs, the threshold is the same as that from the direct filter.

Curve fitting in \code{filterTrended} is done using \code{\link{loessFit}} with a bandwidth of \code{span}.
Lower values may need to be used for a more accurate fit when the trend is highly non-linear.
The bin size is also added to the distance prior to log-transformation, to avoid problems with undefined values when distances are equal to zero.
Empty parts of the interaction space are considered by inferring the abundances and distances of the corresponding bin pairs (though this is skipped if too much of the space is empty).

If \code{reference} is specified, it will be used to compute filter thresholds instead of \code{data}.
This is intended for large bin pairs that have been loaded with \code{filter=1}.
Larger bins provide larger counts for more precise threshold estimates, while the lack of filtering ensures that estimates are not biased.
All threshold estimates are adjusted to account for differences in bin sizes between \code{reference} and \code{data}.
The final values can be used to directly filter on abundances in \code{data}; check out the user's guide for more details.
}

\value{
A list is returned containing \code{abundances}, a numeric vector with the average abundances of all bin pairs in \code{data}.
For \code{filterDirect}, the list contains a numeric scalar \code{threshold}, i.e., the non-specific ligation rate.
For \code{filterTrended}, the list contains \code{threshold}, a numeric vector containing the threshold for each bin pair; and \code{log.distance}, a numeric vector with the log-distances for each bin pair.

If \code{reference} is specified in either function, an additional list named \code{ref} is also returned.
This contains the filtering information for the bin pairs in \code{reference}, same as that reported above for each bin pair in \code{data}.
}

\seealso{
\code{\link{squareCounts}},
\code{\link{scaledAverage}}
}

\author{
Aaron Lun
}

\examples{
# Setting up the object.
a <- 10
b <- 20
regions <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))

set.seed(138153)
npairs <- 500
all.anchors <- sample(length(regions), npairs, replace=TRUE)
all.targets <- as.integer(runif(npairs, 1, all.anchors+1))
y <- DIList(matrix(rnbinom(npairs*4, mu=10, size=10), npairs, 4), anchors=all.anchors, 
    targets=all.targets, regions=regions, exptData=List(width=1))

# Requiring at least 1.5-fold change.
direct <- filterDirect(y)
keep <- direct$abundances > direct$threshold + log2(1.5)
y[keep,]

# Requiring to be above the threshold.
trended <- filterTrended(y)
keep <- trended$abundances > trended$threshold 
y[keep,]

# Running reference comparisons, using larger bin pairs.
w <- 5L
a2 <- a/w
b2 <- b/w
regions2 <- GRanges(rep(c("chrA", "chrB"), c(a2, b2)), 
    IRanges(c(1:a2, 1:b2)*w-w+1L, c(1:a2, 1:b2)*w))
npairs2 <- 20
all.anchors2 <- sample(length(regions2), npairs2, replace=TRUE)
all.targets2 <- as.integer(runif(npairs2, 1, all.anchors2+1))
y2 <- DIList(matrix(rnbinom(npairs2*4, mu=10*w^2, size=10), npairs2, 4), 
    anchors=all.anchors2, targets=all.targets2, regions=regions2, 
    totals=y$totals, exptData=List(width=w))

direct2 <- filterDirect(y, reference=y2)
sum(direct2$abundances > direct2$threshold + log2(1.5))
trended2 <- filterTrended(y, reference=y2)
sum(trended2$abundances > trended2$threshold)
}

\references{
Lin, YC et al. (2012) Global changes in the nuclear positioning of genes and intra- and interdomain genomic interactions that orchestrate B cell fate. \emph{Nat. Immunol.} 13. 1196-1204
}

\keyword{filtering}
