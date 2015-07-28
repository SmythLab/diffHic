\name{boxPairs}
\alias{boxPairs}

\title{Put bin pairs into boxes}
\description{Match smaller bin pairs to the larger bin pairs in which they are nested.}

\usage{
boxPairs(..., reference, minbox=FALSE)
}

\arguments{
\item{...}{one or more named \code{DIList} objects produced by \code{\link{squareCounts}}, with smaller bin sizes than \code{reference}}
\item{reference}{an integer scalar specifying the reference bin size}
\item{minbox}{a logical scalar indicating whether coordinates for the minimum bounding box should be returned}
}

\value{
A list object containing:
\item{indices}{a named list of integer vectors for every \code{DIList} in the ellipsis, see Details} 
\item{anchors,targets}{\code{GRanges} objects specifying the coordinates of the parent bin pair or, if \code{minbox=TRUE}, the minimum bounding box}
}

\details{
Consider the bin size specified in \code{reference}. 
Pairs of these bins are referred to here as the parent bin pairs, and are described in the output \code{pairs} and \code{region}.
The function accepts a number of \code{DIList} objects of bin pair data in the ellipsis, referred to here as input bin pairs.
The aim is to identify the parent bin pair in which each input bin pair is nested.

All input \code{DIList} objects in the ellipsis must be constructed carefully.
In particular, the value of \code{width} in \code{\link{squareCounts}} must be such that \code{reference} is an exact multiple of each \code{width}.
This is necessary to ensure complete nesting.
Otherwise, the behavior of the function will not be clearly defined.

In the output, one vector will be present in \code{indices} for each input \code{DIList} in the ellipsis.
In each vector, each entry represents an index for a single input bin pair in the corresponding \code{DIList}.
This index points to the entries in \code{anchors} and \code{targets} that specify the coordinates of the parent bin pair.
Thus, bin pairs with the same index are nested in the same parent.

Some users may wish to identify bin pairs in one \code{DIList} that are nested within bin pairs in another \code{DIList}.
This can be done by supplying both \code{DIList} objects in the ellipsis, and leaving \code{reference} unspecified.
The value of \code{reference} will be automatically selected as the largest width of the supplied \code{DIList} objects.
Nesting can be identified by \code{\link{match}}ing the output indices for the smaller bin pairs to those of the larger bin pairs.

If \code{minbox=TRUE}, the coordinates in \code{anchors} and \code{targets} represent the minimum bounding box for all nested bin pairs in each parent.
This may be more precise if nesting only occurs in a portion of the interaction space of the parent bin pair.
}

\seealso{
\code{\link{squareCounts}}
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
y <- DIList(matrix(0, ncol(all.combos), 1), anchors=all.combos[2,],
    targets=all.combos[1,], regions=cuts, exptData=List(param=param, width=1))

a5 <- a/5
b5 <- b/5
all.combos2 <- combn(length(cuts)/5, 2) # Bin size of 5.
y2 <- DIList(matrix(0, ncol(all.combos2), 1), anchors=all.combos2[2,],
    targets=all.combos2[1,], exptData=List(param=param, width=5),
	regions=GRanges(rep(c("chrA", "chrB"), c(a5, b5)), 
		IRanges(c((1:a5-1)*5+1, (1:b5-1)*5+1), c(1:a5*5, 1:b5*5))))

# Clustering.
boxPairs(reference=5, larger=y2, smaller=y)
boxPairs(reference=10, larger=y2, smaller=y)
boxPairs(reference=10, larger=y2, smaller=y, minbox=TRUE)
boxPairs(larger=y2, smaller=y)
}

\keyword{clustering}
