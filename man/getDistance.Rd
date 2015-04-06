\name{getDistance}
\alias{getDistance}

\title{Get the linear distance for each interaction}

\description{Compute the distance between interacting regions on the linear genome.}

\usage{
getDistance(data, type=c("mid", "gap", "span"))
}

\arguments{
\item{data}{a \code{DIList} object}
\item{type}{a character string specifying the type of distance to compute}
}

\value{
An integer vector is returned containing the distances between interacting regions for each pair in \code{data}.
}

\details{
The \code{getDistance} function examines each interaction in \code{data$pairs} and computes the distance between the corresponding regions on the linear genome.
Options to compute the distances are \code{mid}, for the distance between the midpoints of the regions; \code{gap}, for the distance of the interval lying between the regions (set to zero, if regions are overlapping); and \code{span}, for the total distance spanned by the interaction, including both regions. 
Interchromosomal interactions are marked with \code{NA}.
}

\author{Aaron Lun}

\seealso{
	\code{\link{squareCounts}},
	\code{\link{connectCounts}}
}

\examples{
# Making up a DIList.
nregs <- 20
reg.sizes <- as.integer(runif(nregs, 5, 10))
ends <- cumsum(reg.sizes)
regions <- GRanges("chrA", IRanges(c(1, ends[-nregs]+1), ends))

all.combos <- combn(length(regions), 2)
y <- DIList(matrix(0, ncol(all.combos), 1), anchors=all.combos[2,],
	targets=all.combos[1,], regions=regions, exptData=List(width=1))

# Collating to count combinations.
getDistance(y)
getDistance(y, type="gap")
getDistance(y, type="span")
}

\keyword{filtering}
