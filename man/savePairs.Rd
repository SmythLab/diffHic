\name{savePairs}
\alias{savePairs}

\title{Save Hi-C interactions}

\description{Save a dataframe of interactions into a directory structure for rapid chromosomal access.}

\usage{
savePairs(x, file, param)
}

\arguments{
	\item{x}{a sorted dataframe with integer fields \code{anchor1.id} and \code{anchor2.id} for each interaction}
	\item{file}{a character string specifying the path for the output index file}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
}

\value{
An index file is produced at the specified \code{file} location, containing the interaction data.
A \code{NULL} value is invisibly returned.
}

\details{
This function facilitates the input of processed Hi-C data from other sources into the current pipeline. 
Each entry in \code{x$anchor1.id} and \code{x$anchor2.id} refers to the index of \code{param$fragments} to denote the interacting regions. 
The \code{x} object will be resorted by \code{anchor1.id}, then \code{anchor2.id}.
If necessary, anchor1 and anchor2 indices will be switched such that the former is never less than the latter. 

The coordinates of the restriction fragment boundaries in \code{param$fragments} should correspond to the reference genome being used. 
In most cases, these can be generated using the \code{\link{cutGenome}} function from any given \code{BSgenome} object.
Values of \code{param$discard} and \code{param$restrict} will not be used here and can be ignored.

Any additional fields in \code{x} will also be saved to file.
Users are recommended to put in \code{anchor1.pos}, \code{anchor1.len}, \code{anchor2.pos} and \code{anchor2.len} fields.
This will allow removal of reads in \code{param$discard} during counting (e.g., with \code{\link{squareCounts}}).
It will also allow proper calculation of statistics with \code{\link{getPairData}}, and quality control with \code{\link{prunePairs}}.
See the output of \code{\link{preparePairs}} for more details on the specification of the fields.

For storing DNase Hi-C data, \code{param$fragments} should be empty, and \code{anchor1.id} and \code{anchor2.id} should contain the indices of the \code{seqlengths}.
Note that both of these fields will be set to zero in the index file - see \code{\link{prepPseudoPairs}} for more details.
In this case, it is essential to store the position and length fields mentioned above.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <-readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts)

n <- 1000
all.a <- as.integer(runif(n, 1L, length(cuts)))
all.t <- as.integer(runif(n, 1L, length(cuts)))
x <- data.frame(anchor1.id=all.a, anchor2.id=all.t,
	anchor1.pos=runif(1:100), anchor1.len=10, 
	anchor2.pos=runif(1:100), anchor2.len=-10)

fout <- "temp2.h5"
savePairs(x, fout, param)
require(rhdf5)
head(h5read(fout, "chrA/chrA"))

\dontshow{
unlink(fout)
}
}

\author{Aaron Lun}

\seealso{
	\code{\link{preparePairs}}, 
	\code{\link{cutGenome}}
}

\keyword{preprocessing}
