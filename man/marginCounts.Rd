\name{marginCounts}
\alias{marginCounts}

\title{Collect marginal counts for each bin}

\description{Count the number of read pairs mapped to each bin across multiple Hi-C libraries.}

\usage{
marginCounts(files, param, width=50000)
}

\arguments{
	\item{files}{a character vector containing paths to the index files}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{width}{an integer scalar specifying the width of each bin} 
}

\value{
A RangedSummarizedExperiment object containing the marginal counts for each bin.
}

\details{
The genome is first split into non-overlapping adjacent bins of size \code{width}, which are rounded to the nearest restriction site. 
The marginal count for each bin is defined as the number of reads in the library mapped to the bin.
This acts as a proxy for genomic coverage by treating Hi-C data as single-end.

Each row of the output RangedSummarizedExperiment refers to a single bin in the linear genome, instead of a bin pair in the interaction space.
The count matrix for all row can be extracted using the \code{assay} method.
Bin coordinates can be extracted using the \code{rowRanges} method.

Larger marginal counts can be collected by increasing the \code{width} value. 
However, this comes at the cost of spatial resolution as adjacent events in the same bin can no longer be distinguished. 
Note that filtering is automatically performed to remove empty bins.

Counting will consider the values of \code{restrict}, \code{discard} and \code{cap} in \code{param}. 
See \code{\link{pairParam}} for more details.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, fout))

# Collating to count combinations.
mar <- marginCounts(fout, param, width=10)
head(assay(mar))
mar <- marginCounts(fout, param, width=50)
head(assay(mar))
mar <- marginCounts(fout, param, width=100)
head(assay(mar))

# Attempting with other parameters.
mar <- marginCounts(fout, reform(param, restrict="chrA"), width=50)
head(assay(mar))
mar <- marginCounts(fout, reform(param, cap=1), width=50)
head(assay(mar))
mar <- marginCounts(fout, reform(param, discard=GRanges("chrA", IRanges(1, 50))), width=50)
head(assay(mar))

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\author{Aaron Lun}

\seealso{
\code{\link{squareCounts}},
\code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
}

\keyword{counting}
