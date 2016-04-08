\name{annotatePairs}
\alias{annotatePairs}

\title{Annotate bin pairs}
\description{Annotate bin pairs based on features overlapping the anchor regions.}

\usage{
annotatePairs(data.list, regions, rnames=names(regions), indices, ...) 
}

\arguments{
\item{data.list}{An InteractionSet or a list of InteractionSet objects containing bin pairs.}
\item{regions}{A GRanges object containing coordinates for the regions of interest.}
\item{rnames}{A character vector containing names to be used for annotation.}
\item{indices}{An integer vector or list of such vectors, indicating the cluster identity for each interaction in \code{data.list}.}
\item{...}{Additional arguments to pass to \code{\link{findOverlaps}}.}
}

\value{
A list of two character vectors \code{anchor1} and \code{anchor2} is returned, containing comma-separated strings of \code{names} for entries in \code{regions} overlapped by the first and second anchor regions respectively.
If \code{indices} is not specified, overlaps are identified to anchor regions of each interaction in \code{data.list}.
Otherwise, overlaps are identified to anchor regions for any interaction in each cluster.
}

\details{
Entries in \code{regions} are identified with any overlap to anchor regions for interactions in \code{data.list}.
The \code{names} for these entries are concatenated into a comma-separated string for easy reporting.
Typically, gene symbols are used in \code{names}, but other values can be supplied depending on the type of annotation.
This is done separately for the first and second anchor regions so that potential interactions between features of interest can be identified.

If \code{indices} is supplied, all interactions corresponding to each unique index are considered to be part of a single cluster.
Overlaps with all interactions in the cluster are subsequently concatenated into a single string.
Cluster indices should range from \code{[1, nclusters]} for any given number of clusters.
This means that the annotation for a cluster corresponding to a certain index can be obtained by subsetting the output vectors with that index.
Otherwise, if \code{indices} is not set, all interactions are assumed to be their own cluster, i.e., annotation is returned for each interaction separately.

Multiple InteractionSet objects can be supplied in \code{data.list}, e.g., if the cluster consists of bin pairs of different sizes.
This means that \code{indices} should also be a list of vectors where each vector indicates the cluster identity of the entries in the corresponding InteractionSet of \code{data.list}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{findOverlaps}},
\code{\link{clusterPairs}}
}

\examples{
# Setting up the objects.
a <- 10
b <- 20
cuts <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))
param <- pairParam(cuts)

all.combos <- combn(length(cuts), 2) 
y <- InteractionSet(matrix(0, ncol(all.combos), 1), 
    GInteractions(anchor1=all.combos[2,], anchor2=all.combos[1,], regions=cuts, mode="reverse"),
    colData=DataFrame(lib.size=1000), metadata=List(param=param, width=1))

regions <- GRanges(rep(c("chrA", "chrB"), c(3,2)), IRanges(c(1,5,8,3,3), c(1,5,8,3,4)))
names(regions) <- LETTERS[seq_along(regions)] 
out <- annotatePairs(y, regions=regions)

# Again, with indices:
indices <- sample(20, length(y), replace=TRUE)
out <- annotatePairs(y, regions=regions, indices=indices)

# Again, with multiple InteractionSet objects:
out <- annotatePairs(list(y, y[1:10,]), regions=regions, indices=list(indices, indices[1:10]))
}

\keyword{annotation}
