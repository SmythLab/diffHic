\name{compartmentalize}
\alias{compartmentalize}

\title{Identify genomic compartments}
\description{Use contact matrices to identify self-interacting genomic compartments}

\usage{
compartmentalize(data, centers=2, dist.correct=TRUE, 
    cov.correct=TRUE, robust.cov=5, inter=FALSE, ...)
}

\arguments{
\item{data}{a \code{DIList} object containing bin pair data, like that produced by \code{\link{squareCounts}}}
\item{centers}{an integer scalar, specifying the number of clusters to form in \code{\link{kmeans}}}
\item{dist.correct}{a logical scalar, indicating whether abundances should be corrected for distance biases}
\item{cov.correct}{a logical scalar, indicating whether abundances should be corrected for coverage biases}
\item{robust.cov}{a numeric scalar, specifying the multiple of MADs beyond which coverage outliers are removed}
\item{inter}{a logical scalar, indicating whether inter-chromosomal space should be used during clustering}
\item{...}{other arguments to pass to \code{\link{kmeans}}}
}

\details{
This function uses the interaction space to partition each linear chromosome into compartments.
Bins in the same compartment interact more frequently with each other compared to bins in different compartments.
This forms a checkerboard-like pattern in the interaction space that can be used to define the genomic intervals in each compartment.
Typically, one compartment is gene-rich and is defined as ``open'', while the other is gene-poor and defined as ``closed''.

Compartment identification is done by setting up a contact matrix where each row/rolumn represents a bin and each matrix entry contains the frequency of contacts between bins.
Bins (i.e., rows) with similar interaction profiles (i.e., entries across columns) are clustered together with the k-means method.
Those with the same ID in the output \code{compartment} vector are in the same compartment.

If \code{dist.correct=TRUE}, frequencies are normalized to mitigate the effect of distance and to improve the visibility of long-range interactions.
This is done by computing the residuals of the distance-dependent trend - see \code{\link{filterTrended}} for more details.
If \code{cov.correct=TRUE}, frequencies are also normalized to eliminate coverage biases betwen bins.
This is done by computing the average coverage of each row/column, and dividing each matrix entry by the square root averages of the relevant row and column.

Extremely low-coverage regions such as telomeres and centromeres can confound k-means clustering.
To protect against this, all bins with (distance-corrected) coverages that are more than \code{robust.cov} MAD's away from the median coverage of each chromosome are identified and removed.
These bins will be marked with \code{NA} in the returned \code{compartment} for that chromosome.
To turn off robustification, set \code{robust.cov} to \code{NA}.

When \code{inter=FALSE}, clustering is performed for each chromosome separately, using a series of intra-chromosomal contact matrices.
This saves memory and performs well, given that most information is contained in the intra-chromosomal space.
If \code{inter} is set to \code{TRUE}, clustering is performed for the entire genome, including the interaction space between different chromosomes, i.e., the inter-chromomsal space.
This may be useful for matching compartments of the same type between chromosomes, but will use a lot more memory.

By default, \code{centers} is set to 2 to model the open and closed compartments.
While a larger value can be used to obtain more clusters, care is required as the interpretation of the resulting compartments becomes more difficult.
If desired, users can also apply their own clustering methods on the \code{matrix} returned in the output.
}

\value{
If \code{inter=FALSE}, a named list of lists is returned where each internal list corresponds to a chromosome in \code{data} and contains \code{compartment}, an integer vector of compartment IDs for all bins in that chromosome; and \code{matrix}, a numeric matrix of (normalized) contact frequencies for the intra-chromosomal space. 
Entries in \code{compartment} and rows/columns in \code{matrix} are named according to the matching index of \code{regions(data)}.

If \code{inter=TRUE}, a list is returned that contains \code{compartment}, an integer vector of IDs for all bins in the genome; and \code{matrix}, a contact matrix spanning the entire interaction space.
Each entry of \code{compartment} and each row/column in \code{matrix} corresponds to an entry of \code{regions(data)}.
}

\seealso{
\code{\link{squareCounts}},
\code{\link{filterTrended}},
\code{\link{kmeans}}
}

\author{Aaron Lun}

\examples{
# Dummying up some data.
set.seed(3426)
npts <- 100
npairs <- 5000
nlibs <- 4
anchors <- sample(npts, npairs, replace=TRUE)
targets <- sample(npts, npairs, replace=TRUE)
data <- DIList(counts=matrix(rpois(npairs*nlibs, runif(npairs, 10, 100)), nrow=npairs),
   totals=runif(nlibs, 1e6, 2e6), anchors=pmax(anchors, targets), targets=pmin(anchors, targets),
   regions=GRanges(c(rep("chrA", 80), rep("chrB", 20)), IRanges(c(1:80, 1:20), c(1:80, 1:20))),
   exptData=List(width=1))

# Running compartmentalization.
out <- compartmentalize(data)
head(out$chrA$compartment)
dim(out$chrA$matrix)
head(out$chrB$compartment)
dim(out$chrB$matrix)

out2 <- compartmentalize(data, inter=TRUE)
head(out2$compartment)
out2$matrix[1,]

test <- compartmentalize(data, cov.correct=FALSE)
test <- compartmentalize(data, dist.correct=FALSE)
test <- compartmentalize(data, robust.cov=NA)
}

\keyword{clustering}

\references{
Lieberman-Aiden E et al. (2009). Comprehensive mapping of long-range interactions reveals folding principles of the human genome. \emph{Science} 326, 289-293. 

Lajoie BR, Dekker J, Kaplan N (2014). The hitchhiker's guide to Hi-C analysis: practical guidelines. \emph{Methods} 72, 65-75.
}
