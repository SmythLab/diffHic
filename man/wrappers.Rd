\name{DIList-wrappers}
\alias{asDGEList}
\alias{asDGEList,DIList-method}
\alias{normalize}
\alias{normalize,DIList-method}

\title{Statistical wrappers for DIList objects}

\description{
Convenience wrappers for statistical routines operating on DIList objects.
}

\usage{
\S4method{normalize}{DIList}(object, lib.sizes, ...)
\S4method{asDGEList}{DIList}(object, lib.sizes, ...)
}

\arguments{
	\item{object}{a \code{DIList} object, like that produced by \code{\link{squareCounts}}}
	\item{lib.sizes}{an (optional) integer vector of library sizes}
	\item{...}{other arguments to be passed to the function being wrapped}
}

\details{
Counts are extracted from the matrix in the \code{DIList} object.
If not specified in \code{lib.sizes}, library sizes are taken from the \code{totals} field in the column data of \code{object}.
Warnings will be generated if this field is not present.

In the \code{normalize} method, the extracted counts and library sizes are supplied to \code{\link{normalizeCounts}}.
Similarly, the \code{asDGEList} method wraps the \code{\link{DGEList}} constructor.
In both cases, any arguments in \code{...} are also passed to the wrapped functions.
}

\value{
For \code{normalize}, either a numeric matrix or vector is returned; see \code{\link{normalizeCounts}}.

For \code{asDGEList}, a \code{DGEList} object is returned.
}

\author{Aaron Lun}

\examples{
blah <- DIList(counts=matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), ncol=2), 
    totals=c(10L, 10L), anchors=c(1,2,3,4,5,6), targets=c(1,1,2,2,3,3), 
	regions=GRanges("chrA", IRanges(10+1:20, 2+21:40)))

asDGEList(blah)
asDGEList(blah, lib.sizes=c(1,2)) # Manually set to some other value.
asDGEList(blah, group=c("a", "b")) 

normalize(blah)  
normalize(blah, lib.sizes=c(1,2)) # Manually set to some other value. 
normalize(blah, logratioTrim=0)
normalize(blah, type="loess")

c(blah, blah)
c(blah[1:2,], blah[3:4,], add.totals=FALSE)
}

\seealso{
\code{\link{DGEList}},
\code{\link{normalizeCounts}},
\code{\link{squareCounts}}
}

\keyword{normalization}
