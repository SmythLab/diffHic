\name{boxPairs}
\alias{boxPairs}

\title{Put bin pairs into boxes}
\description{Match smaller bin pairs to the larger bin pairs in which they are nested.}

\usage{
boxPairs(..., reference, minbox=FALSE, index.only=FALSE)
}

\arguments{
\item{...}{One or more named InteractionSet objects produced by \code{\link{squareCounts}}, with smaller bin sizes than \code{reference}.}
\item{reference}{An integer scalar specifying the reference bin size.}
\item{minbox}{A logical scalar indicating whether coordinates for the minimum bounding box should be returned.}
\item{index.only}{A logical scalar indicating whether only indices should be returned.}
}

\value{
If \code{index.only=FALSE}, a named list is returned containing:
\describe{
\item{\code{indices}:}{a named list of integer vectors for every InteractionSet in the ellipsis, see Details.} 
\item{\code{interactions}:}{A ReverseStrictGInteractions object containing the coordinates of the parent bin pair or, if \code{minbox=TRUE}, the minimum bounding box.}
}
If \code{index.only=TRUE}, the \code{indices} are returned directly without computing coordinates.
}

\details{
Consider the bin size specified in \code{reference}. 
Pairs of these bins are referred to here as the parent bin pairs, and are described in the output \code{pairs} and \code{region}.
The function accepts a number of InteractionSet objects of bin pair data in the ellipsis, referred to here as input bin pairs.
The aim is to identify the parent bin pair in which each input bin pair is nested.

All input InteractionSet objects in the ellipsis must be constructed carefully.
In particular, the value of \code{width} in \code{\link{squareCounts}} must be such that \code{reference} is an exact multiple of each \code{width}.
This is necessary to ensure complete nesting.
Otherwise, the behavior of the function will not be clearly defined.

In the output, one vector will be present in \code{indices} for each input InteractionSet in the ellipsis.
In each vector, each entry represents an index for a single input bin pair in the corresponding InteractionSet.
This index points to the entries in \code{interactions} that specify the coordinates of the parent bin pair.
Thus, bin pairs with the same index are nested in the same parent.

Some users may wish to identify bin pairs in one InteractionSet that are nested within bin pairs in another InteractionSet.
This can be done by supplying both InteractionSet objects in the ellipsis, and leaving \code{reference} unspecified.
The value of \code{reference} will be automatically selected as the largest width of the supplied InteractionSet objects.
Nesting can be identified by \code{\link{match}}ing the output indices for the smaller bin pairs to those of the larger bin pairs.

If \code{minbox=TRUE}, the coordinates in \code{interactions} represent the minimum bounding box for all nested bin pairs in each parent.
This may be more precise if nesting only occurs in a portion of the interaction space of the parent bin pair.

If \code{index.only=TRUE}, only the indices are returned and coordinates are not computed. 
This is largely for efficiency purposes when \code{boxPairs} is called by internal functions.
}

\seealso{
\code{\link{squareCounts}},
\code{\link{clusterPairs}}
}

\author{
Aaron Lun
}

\examples{
# Setting up the objects.
a <- 10
b <- 20
cuts <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))
param <- pairParam(cuts)

all.combos <- combn(length(cuts), 2) # Bin size of 1.
y <- InteractionSet(matrix(0, ncol(all.combos), 1), 
    GInteractions(anchor1=all.combos[2,], anchor2=all.combos[1,], regions=cuts, mode="reverse"),
    colData=DataFrame(lib.size=1000), metadata=List(param=param, width=1))

a5 <- a/5
b5 <- b/5
all.combos2 <- combn(length(cuts)/5, 2) # Bin size of 5.
y2 <- InteractionSet(matrix(0, ncol(all.combos2), 1), 
    GInteractions(anchor1=all.combos2[2,], anchor2=all.combos2[1,], 
    	regions=GRanges(rep(c("chrA", "chrB"), c(a5, b5)), 
    		IRanges(c((1:a5-1)*5+1, (1:b5-1)*5+1), c(1:a5*5, 1:b5*5))), mode="reverse"),
    colData=DataFrame(lib.size=1000), metadata=List(param=param, width=5))

# Clustering.
boxPairs(reference=5, larger=y2, smaller=y)
boxPairs(reference=10, larger=y2, smaller=y)
boxPairs(reference=10, larger=y2, smaller=y, minbox=TRUE)
boxPairs(larger=y2, smaller=y)
}

\keyword{clustering}
