\name{filterPeaks}
\alias{filterPeaks}

\title{Filter bin pairs for likely peaks}
\description{Identify bin pairs that are likely to represent punctate peaks in the interaction space.}

\usage{
filterPeaks(data, enrichment, assay.bp=1, assay.neighbors=NULL, get.enrich=FALSE,
    min.enrich=log2(1.5), min.count=5, min.diag=2L, ...)
}

\arguments{
\item{data}{an InteractionSet object produced by \code{\link{enrichedPairs}} or \code{\link{neighborCounts}}}
\item{enrichment}{a numeric vector of enrichment values}  
\item{assay.bp}{a string or integer scalar specifying the assay containing bin pair counts}
\item{assay.neighbors}{a character vector containing names for the neighborhood regions, see \code{\link{enrichedPairs}} for details}
\item{get.enrich}{a logical scalar indicating whether enrichment values should be returned}
\item{min.enrich}{a numeric scalar indicating the minimum enrichment score for a peak}
\item{min.count}{a numeric scalar indicating the minimum average count for a peak}
\item{min.diag}{an integer scalar specifying the minimum diagonal in the interaction space with which to consider a peak}
\item{...}{other arguments to be passed to \code{\link{aveLogCPM}} for the average count filter}
}

\value{
If \code{get.enrich=TRUE}, a numeric vector of enrichment values for each bin pair.
Otherwise, a logical vector indicating whether or not each bin pair is to be considered as a peak.
}

\details{
Filtering on the local enrichment scores identifies high-intensity islands in the interaction space.
However, this alone is not sufficient to identify sensible peaks.
Filtering on the absolute average counts prevents the calling of low-abundance bin pairs with high enrichment scores due to empty neighborhoods.
Filtering on the diagonals prevents calling of high-abundance short-range interactions that are usually uninteresting.
If either \code{min.count} or \code{min.diag} are \code{NULL}, no filtering will be performed on the average counts and diagonals, respectively.

To compute enrichment values, we assume that the number of read pairs in neighborhood areas have been counted using \code{\link{enrichedPairs}} or \code{\link{neighborCounts}}.
For a given bin pair in \code{data}, this function computes the mean abundance across libraries for each surrounding neighborhood, scaled by the neighborhood area (i.e., the number of bin pairs it contains).
The local background for the target bin pair is defined as the maximum of the mean abundances for all neighborhoods.
The enrichment value is then defined as the the difference between the target bin pair's abundance and its local background.
The idea is that bin pairs with high enrichments are likely to represent punctate interactions between clearly defined loci.
Selecting for high enrichments can then select for these peak-like features in the interaction space.

% If a peak is outside a TAD (but still within flank distance), then the TAD will be partially included in the 
% background. This could inflate the background abundance and prevent calling of the peak. To avoid this, you'd
% have to be able to specify the boundaries of the structures in the interaction space, which is complicated.
% However, most peaks seem to occur at the cap of domains anyway (or, most looping interactions generate their
% own domain) so it might not be likely to observe a peak outside of a domain. Ah well.

The maximizing strategy is designed to mitigate the effects of structural features.
Region B will capture the high interaction intensity within genomic domains like TADs, while the C and D will capture any bands in the interaction space.
The abundance will be high for any neighborhood that captures a high-intensity feature, as the average counts will be large for all bin pairs within the features.
This will then be chosen as the maximum during calculation of enrichment values.
Otherwise, if only region A were used, the background abundance would be decreased by low-intensity bin pairs outside of the features.
This results in spuriously high enrichment values for target bin pairs on the feature boundaries.

% By default, nothing is done to adjust for the effect of distance on abundance for intra-chromosomal bin pairs.
% This is because the counts are generally too low to routinely fit a reliable trend.
% That said, users can still supply distance-adjusted abundances as \code{abundances}.
% Such values can be defined as the residuals of the fit from \code{\link{filterTrended}}.
% Obviously, no such work is required for inter-chromosomal bin pairs. % Distance adjustment will have no effect, as the direct filtering threshold is constant for all bin pairs.
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
data <- InteractionSet(matrix(as.integer(rnbinom(200, mu=10, size=10)), 50, 4), 
    GInteractions(anchor1=all.anchor1, anchor2=all.anchor2, regions=regions, mode="reverse"), 
    colData=DataFrame(totals=runif(4, 1e6, 2e6)), metadata=List(width=1))

# Getting peaks.
enrichment <- enrichedPairs(data)
summary(filterPeaks(enrichment, min.enrich=0.5))
summary(filterPeaks(enrichment, min.enrich=0.5, min.count=10))
summary(filterPeaks(enrichment, min.enrich=0.5, min.diag=NULL))
}
