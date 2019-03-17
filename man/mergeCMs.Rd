\alias{mergeCMs}
\name{mergeCMs}

\title{Merge ContactMatrix objects}
\description{Merge ContactMatrix objects into an InteractionSet object containing counts for pairs of interacting regions.}

\usage{
mergeCMs(..., deflate.args=list())
}

\arguments{
\item{...}{ContactMatrix objects containing read pair counts for the same area of the interaction space, usually defined as pairs of bins along a genomic interval.}
\item{deflate.args}{A list of arguments to pass to \code{\link[InteractionSet]{deflate}}.}
}

\details{
This function facilitates the conversion of multiple ContactMatrix objects into a single InteractionSet object.
Each ContactMatrix corresponds to a sample and should contain read pair counts between the corresponding row/column regions.
The dimensions of all ContactMatrix objects should be the same and the rows/columns should represent the same genomic regions.

The idea is to produce an object equivalent to the output of \code{\link{squareCounts}} when contact matrices are available instead of BAM files.
This is done via the \code{\link[InteractionSet]{deflate}} method, where each ContactMatrix is converted to an InteractionSet using \code{deflate.args}.
Entries of the ContactMatrix are equivalent to paired regions in an InteractionSet (which, in most cases, are bins of constant width).

The InteractionSet objects for all supplied samples are then combined into a single object for downstream input.
This step will throw errors if the original ContactMatrix objects do not cover the same area of the interaction space.
Column names are set to any names for \code{...}, and the total number of read pairs in each ContactMatrix is stored in \code{totals}.

The \code{width} value in the \code{metadata} of the output InteractionSet is set to the median width of the interacting regions.
The \code{totals} field in the output \code{colData} is also set to be equal to the sum of the counts in each ContactMatrix (after removing redundant regions).
Note that this only makes sense if the ContactMatrix objects contain interactions between non-overlapping genomic bins.

The \code{param} value is not set in the \code{metadata} of the output object.
This depends on how the ContactMatrix objects were constructed in the first place, which is not known to the function.
}

\value{
An InteractionSet object containing counts for interacting regions.
}

\author{
Aaron Lun
}

\seealso{
    \code{\link[InteractionSet]{deflate}},
    \code{\link{squareCounts}},
    \code{\link{connectCounts}}
} 

\examples{
example(ContactMatrix, echo=FALSE)
mergeCMs(x, x2)
}

