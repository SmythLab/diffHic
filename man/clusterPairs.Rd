\name{clusterPairs}
\alias{clusterPairs}

\title{Cluster bin pairs}
\description{Aggregate bin pairs into local clusters for summarization.}

\usage{
clusterPairs(..., tol, upper=1e6, index.only=FALSE)
}

\arguments{
\item{...}{One or more InteractionSet objects, optionally named.}
\item{tol}{A numeric scalar specifying the maximum distance between bin pairs in base pairs.}
\item{upper}{A numeric scalar specifying the maximum size of each cluster in base pairs.}
\item{index.only}{A logical scalar indicating whether only indices should be returned.}

}

\value{
If \code{index.only=FALSE}, a named list is returned containing:
\describe{
\item{\code{indices}:}{A named list of integer vectors where each vector contains a cluster ID for each interaction in the corresponding input InteractionSet object.}
\item{\code{interactions}:}{A ReverseStrictGInteractions object containing the coordinates of the sides of the bounding box for each cluster.}
}
If \code{index.only=TRUE}, the \code{indices} are returned directly without computing coordinates.
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
The coordinates of the box for each cluster is stored in each row of the output \code{interactions}.
The cluster ID in each \code{indices} vector represents the row index for these coordinates.

If \code{index.only=TRUE}, only the indices are returned and coordinates of the bounding box are not computed. 
This is largely for efficiency purposes when \code{clusterPairs} is called by internal functions.
}

\seealso{
\code{\link{squareCounts}},
\code{\link{diClusters}},
\code{\link{boxPairs}}
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
all.anchor1 <- sample(length(regions), 50, replace=TRUE)
all.anchor2 <- as.integer(runif(50, 1, all.anchor1+1))
y <- InteractionSet(matrix(0, 50, 1), 
    GInteractions(anchor1=all.anchor1, anchor2=all.anchor2, regions=regions, mode="reverse"),
    colData=DataFrame(lib.size=1000), metadata=List(width=1))

# Clustering; note, small tolerances are used in this toy example.
clusterPairs(y, tol=1)
clusterPairs(y, tol=3)
clusterPairs(y, tol=5)
clusterPairs(y, tol=5, upper=5)

# Multiple bin sizes allowed.
a2 <- a/2
b2 <- b/2
rep.regions <- GRanges(rep(c("chrA", "chrB"), c(a2, b2)), 
    IRanges(c(1:a2*2, 1:b2*2), c(1:a2*2, 1:b2*2)))
rep.anchor1 <- sample(length(rep.regions), 10, replace=TRUE)
rep.anchor2 <- as.integer(runif(10, 1, rep.anchor1+1))
y2 <- InteractionSet(matrix(0, 10, 1), 
    GInteractions(anchor1=rep.anchor1, anchor2=rep.anchor2, regions=rep.regions, mode="reverse"), 
    colData=DataFrame(lib.size=1000), metadata=List(width=2))

clusterPairs(y, y2, tol=1)
clusterPairs(y, y2, tol=3)
clusterPairs(y, y2, tol=5)
clusterPairs(y, tol=5, upper=5)
}

\keyword{clustering}
