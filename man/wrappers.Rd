\name{DIList-wrappers}
\alias{asDGEList}
\alias{asDGEList,DIList-method}
\alias{normOffsets}
\alias{normOffsets,DIList-method}
\alias{normalize}
\alias{normalize,DIList-method}

\title{Statistical wrappers for DIList objects}

\description{
Convenience wrappers for statistical routines operating on DIList objects.

\bold{NOTE:} These methods are all deprecated, along with the DIList class itself.
Users can apply the same methods to InteractionSet objects.
}

\usage{
\S4method{normOffsets}{DIList}(object, lib.sizes, ...)
\S4method{normalize}{DIList}(object, lib.sizes, ...) # deprecated, use normOffsets
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

In the \code{normOffsets} method, the extracted counts and library sizes are supplied to \code{\link{normOffsets,matrix-method}}.
Similarly, the \code{asDGEList} method wraps the \code{\link{DGEList}} constructor.
In both cases, any arguments in \code{...} are also passed to the wrapped functions.
}

\value{
For \code{normOffsets}, either a numeric matrix or vector is returned; see \code{\link{normOffsets,matrix-method}}.

For \code{asDGEList}, a \code{DGEList} object is returned.
}

\author{Aaron Lun}

