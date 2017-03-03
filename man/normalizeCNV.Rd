\name{normalizeCNV}
\alias{normalizeCNV}
\alias{matchMargins}

\title{Normalize CNV biases}
\description{Compute normalization offsets to remove CNV-driven and abundance-dependent biases}

\usage{
normalizeCNV(data, margins, prior.count=3, span=0.3, maxk=500, 
    assay.data=1, assay.marg=1, ...)

matchMargins(data, margins)
}

\arguments{
\item{data}{an InteractionSet object produced by \code{\link{squareCounts}}}
\item{margins}{a RangedSummarizedExperiment object produced by \code{\link{marginCounts}}}
\item{prior.count}{a numeric scalar specifying the prior count to use in computing marginal log-ratios}
\item{span}{a numeric scalar between 0 and 1, describing the span of the fit}
\item{maxk}{a integer scalar specifying the number of vertices to use during local fitting}
\item{assay.data}{a string or integer scalar specifying the matrix to use from \code{data}}
\item{assay.marg}{a string or integer scalar specifying the matrix to use from \code{margins}}
\item{...}{other arguments to pass to \code{\link{locfit}}}
}

\details{
Each bin pair in \code{data} is associated with three covariates.
The first two are the marginal log-ratios of the corresponding bins, i.e., the log-ratio of the marginal counts between two libraries.
These represent the relative CNVs in the interacting regions between libraries.
To avoid redundancy, the first covariate is the larger marginal log-ratio whereas the second is the smaller.
The third covariate is the average abundance across all libraries.

Each bin pair is also associated with a response, i.e., the log-ratio of the interaction counts between two libraries.
A loess-like surface is fitted to the response against the three covariates, using the \code{\link{locfit}} function.
The aim is to eliminate systematic differences between libraries at any combination of covariate values.
This removes CNV-driven biases as well as trended biases with respect to the abundance.
The fitted value can then be used as a GLM offset for each bin pair.

% CNVs in both interacting regions are used as separate covariates, rather than using the product of the relative copies from the two regions as a single covariate.
% This is because a doubling in copy number for the two regions mightn't have the same effect on intensity as a 4-fold increase in one region. 
% The effect of CNVs will also depend on the mechanism through which ligation events are generated.
% For example, a doubling in copy number for both regions might result in a 4-fold increase for random ligation, but only a 2-fold increase for genuine interactions (as the interaction is just doubled).
% Using the abundance as a covariate will ensure that these cases are fitted separately.

The objects in \code{data} and \code{margins} should be constructed with the same \code{width} and \code{param} arguments in their respective functions.
This ensures that the regions are the same, so that the marginal counts can be directly used.
Matching of the bins in each bin pair in \code{data} to indices of \code{margins} is performed using \code{matchMargins}.
Note that the marginal counts are not directly computed from \code{data} as filtering of bin pairs may be performed beforehand.

In practice, normalization offsets are computed for each library relative to a single reference ``average'' library.
This average library is constructed by using the average abundance as the (log-)count for both the bin pair and marginal counts.
The space of all pairs of CNV log-ratios is also rotated by 45 degrees prior to smoothing.
This improves the performance of the approximations used by \code{\link{locfit}}.

The fit parameters can be changed by varying \code{span}, \code{maxk} and additional arguments in \code{\link{locfit}}.
Higher values of \code{span} will increase smoothness, at the cost of sensitivity.
Increases in \code{maxk} may be required to obtain a more accurate approximation when fitting large datasets.
In all cases, a loess fit of degree 1 is used.

For use by downstream functions, the offset matrix can be stored as an assay named \code{"offset"} in the \code{data} object.
}

\value{
For \code{normalizeCNV}, a numeric matrix is returned with the same dimensions as \code{counts(data)}.
This contains log-based GLM offsets for each bin pair in each library.

For \code{matchMargins}, a data frame is returned with integer fields \code{anchor1} and \code{anchor2}.
Each field specifies the index in \code{margins} corresponding to the bin regions for each bin pair in \code{data}. 
}

\seealso{
\code{\link{locfit}}, 
\code{\link{lp}},
\code{\link{squareCounts}},
\code{\link{marginCounts}}
}

\author{
Aaron Lun
}

\examples{
# Dummying up some data.
set.seed(3423746)
npts <- 100
npairs <- 5000
nlibs <- 4
anchor1 <- sample(npts, npairs, replace=TRUE)
anchor2 <- sample(npts, npairs, replace=TRUE)

data <- InteractionSet(
    list(counts=matrix(rpois(npairs*nlibs, runif(npairs, 10, 100)), nrow=npairs)),
    GInteractions(anchor1=anchor1, anchor2=anchor2, 
        regions=GRanges("chrA", IRanges(1:npts, 1:npts)), mode="reverse"),
	colData=DataFrame(totals=runif(nlibs, 1e6, 2e6))) 

margins <- SummarizedExperiment(matrix(rpois(npts*nlibs, 100), nrow=npts),
	colData=DataFrame(totals=data$totals), rowRanges=regions(data))

# Running normalizeCNV.
out <- normalizeCNV(data, margins)
head(out)
head(normalizeCNV(data, margins, prior.count=1))
head(normalizeCNV(data, margins, span=0.5))

# Store offsets as the 'offset' assay for use by, e.g., asDGEList.
assays(data)$offset <- out
data

# Occasionally locfit will complain; increase maxk to compensate.
data <- InteractionSet(matrix(rpois(npairs*nlibs, 20), nrow=npairs),
    GInteractions(anchor1=anchor1, anchor2=anchor2,
        regions=GRanges("chrA", IRanges(1:npts, 1:npts)), mode="reverse"),
	colData=DataFrame(totals=runif(nlibs, 1e6, 2e6)))
tryCatch(head(normalizeCNV(data, margins, maxk=100)), error=function(e) e)
head(normalizeCNV(data, margins, maxk=1000))

# Matching margins.
matched <- matchMargins(data, margins)
head(matched)
anchor1.counts <- margins[matched$anchor1,]
anchor2.counts <- margins[matched$anchor2,]
}

\keyword{normalization}
