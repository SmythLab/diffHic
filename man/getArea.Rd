\name{getArea}
\alias{getArea}

\title{Get interaction area}

\description{Compute area in the interaction space for each pair of regions.}

\usage{
getArea(data, bp=TRUE)
}

\arguments{
\item{data}{a \code{DIList} object}
\item{bp}{a logical scalar indicating whether areas should be reported in base-pair terms}
}

\value{
A numeric vector is returned containing the area in the interaction space for each pair of regions in \code{data}.
}

\details{
The \code{getArea} function returns the area in the interaction space for each pair of regions.
If \code{bp=TRUE}, the area is reported in terms of squared base pairs.
This tends to be the easiest to interpret.
Otherwise, the area is reported as the number of pairs of restriction fragments.
This may be more relevant to the actual resolution of the Hi-C experiment.

Some special consideration is required for areas overlapping the diagonal.
This is because counting is only performed on one side of the diagonal, to avoid redundancy.
Base-pair areas are automatically adjusted to account for this feature, based on the presence of partial overlaps between interacting regions.

For fragment-based areas, some additional work is required to properly compute areas around the diagonal for partially overlapping regions.
This is only necessary when \code{data} is produced by \code{\link{connectCounts}}.
This is because bins will not partially overlap in any significant manner when counts are generated with \code{\link{squareCounts}}.
}

\author{Aaron Lun}

\seealso{
	\code{\link{squareCounts}},
	\code{\link{connectCounts}}
}

\examples{
# Making up a DIList for binned data.
nfrags <- 50
frag.sizes <- as.integer(runif(nfrags, 5, 10))
ends <- cumsum(frag.sizes)
cuts <- GRanges("chrA", IRanges(c(1, ends[-nfrags]+1), ends))
param <- pairParam(cuts) 

regions <- diffHic:::.getBinID(cuts, 20)$region
all.combos <- combn(length(regions), 2)
y <- DIList(matrix(0, ncol(all.combos), 1), anchors=all.combos[2,],
    targets=all.combos[1,], regions=regions, exptData=List(param=param, width=20))

# Generating partially overlapping regions.
set.seed(3424)
re <- sample(nfrags, 20)
rs <- as.integer(runif(20, 1, re+1))
regions <- GRanges("chrA", IRanges(start(cuts)[rs], end(cuts)[re]))
regions$nfrags <- re - rs + 1L
regions <- sort(regions)
all.combos <- combn(length(regions), 2)
y2 <- DIList(matrix(0, ncol(all.combos), 1), anchors=all.combos[2,],
    targets=all.combos[1,], regions=regions, exptData=List(param=param))

#### Getting areas. ####

getArea(y)
getArea(y, bp=FALSE)

getArea(y2)
getArea(y2, bp=FALSE)
}

\keyword{filtering}
