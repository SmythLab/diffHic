\name{clusterPairs}
\alias{clusterPairs}

\title{Cluster bin pairs}
\description{Aggregate bin pairs into local clusters for summarization.}

\usage{
clusterPairs(..., tol, upper=1e6)
}

\arguments{
\item{...}{one or more \code{DIList} objects, optionally named}
\item{tol}{a numeric scalar specifying the maximum distance between bin pairs}
\item{upper}{a numeric scalar specifying the maximum size of each cluster}
}

\value{
A list containing \code{indices}, a named list of integer vectors where each vector contains a cluster ID for each interaction in the corresponding input \code{DIList} object;
\code{anchors}, a \code{GRanges} object containing the dimensions of the bounding box of each cluster on the anchor chromosome;
and \code{targets}, a \code{GRanges} object containing the corresponding dimensions on the target chromosome;

}

\details{
Clustering is performed by putting a interaction in a cluster if the smallest Chebyshev distance to any interaction already inside the cluster is less than \code{tol}. 
This is a cross between single-linkage approaches and density-based methods, especially after filtering removes low-density regions.
In this manner, adjacent events in the interaction space can be clustered together.
Interactions that are assigned with the same cluster ID belong to the same cluster.

The input \code{data} objects can be taken from the output of \code{\link{squareCounts}} or \code{\link{connectCounts}}.
For the former, inputs can consist of interactions with multiple bin sizes.
It would be prudent to filter the former based on the average abundances, to reduce the density of bin pairs in the interaction space. 
Otherwise, clusters may be too large to be easily interpreted.

Alternatively, to avoid excessively large clusters, this function can also split each cluster into roughly-equally sized subclusters.
The maximum value of any dimension of the subclusters is approxiamtely equal to \code{upper}. 
This aims to improve the spatial interpretability of the clustering result.

% Having evenly sized clusters provides best spatial resolution and interpretability, e.g., much like binary searching.
% However, it's not clear how you can provide an even size while enforcing a strict upper bound of `upper`.

There is no guarantee that each cluster forms a regular shape in the interaction space.
Instead, a minimum bounding box is reported containing all bin pairs in each cluster.
The coordinates of the box for each cluster is stored in each row of the output \code{anchors} and \code{targets}.
The cluster ID in each \code{indices} vector represents the row index for these coordinates.
}

\seealso{
\code{\link{squareCounts}}
}

\author{
Aaron Lun
}

\examples{
# Setting up the object.
a <- 10
b <- 20
regions <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))

set.seed(3423)
all.anchors <- sample(length(regions), 50, replace=TRUE)
all.targets <- as.integer(runif(50, 1, all.anchors+1))
y <- DIList(matrix(0, 50, 1), anchors=all.anchors, targets=all.targets, 
    regions=regions, exptData=List(width=1))

# Clustering; note, small tolerances are used in this toy example.
clusterPairs(y, tol=1)
clusterPairs(y, tol=3)
clusterPairs(y, tol=5)
clusterPairs(y, tol=5, upper=5)

# Multiple bin sizes allowed.
a2 <- a/2
b2 <- b/2
regions2 <- GRanges(rep(c("chrA", "chrB"), c(a2, b2)), 
    IRanges(c(1:a2*2, 1:b2*2), c(1:a2*2, 1:b2*2)))
all.anchors2 <- sample(length(regions2), 10, replace=TRUE)
all.targets2 <- as.integer(runif(10, 1, all.anchors2+1))
y2 <- DIList(matrix(0, 10, 1), anchors=all.anchors2, targets=all.targets2, 
    regions=regions2, exptData=List(width=2))

clusterPairs(y, y2, tol=1)
clusterPairs(y, y2, tol=3)
clusterPairs(y, y2, tol=5)
clusterPairs(y, tol=5, upper=5)
}

\keyword{clustering}
