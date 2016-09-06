\name{domainDirections}
\alias{domainDirections}

\title{Calculate domain directionality}
\description{Collect directionality statistics for domain identification with genomic bins.}

\usage{
domainDirections(files, param, width=50000, span=10)
}

\arguments{
\item{files}{a character vector containing paths to the index files generated from each Hi-C library}
\item{param}{a \code{pairParam} object containing read extraction parameters}
\item{width}{an integer scalar specifying the width of each bin in base pairs}
\item{span}{an integer scalar specifying the distance to consider for up/downstream interactions}
}

\details{
The genome is partitioned into bins of size \code{width}.
For each bin, this function computes the total number of read pairs between that bin and the \code{span} upstream bins (i.e., those with higher genomic coordinates).
This is repeated for the \code{span} downstream bins, thus yielding two counts (up and down) per bin.

A RangedSummarizedExperiment is returned containing the coordinates of each bin and two matrices of counts, named \code{"up"} and \code{"down"}.
Each row of the matrix corresponds to a bin, while each column corresponds to a library in \code{files}.
Each entry of the matrix stores the total count of read pairs to upstream or downstream bins.

The total up- and downstream counts can be used to compute a directionality statistic, e.g., as defined by Dixon et al, or by computing the log-fold change between fields.
Alternatively, it can be used to identify differential domains - see the user's guide for more details.
}

\value{
A RangedSummarizedExperiment object with one row for each bin in the genome.
It contains two integer matrices named \code{"up"} and \code{"down"}, containing the counts to upstream and downstream bins respectively.
}

\seealso{
\code{\link{squareCounts}}
}

\author{Aaron Lun}

\references{
Dixon JR et al. (2012). Topological domains in mammalian genomes identified by analysis of chromatin interactions. \emph{Nature} 485:376-380.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, file=fout))

# Not really that informative; see user's guide.
out <- domainDirections(fout, param, width=10)
out

# Calculating directionality log-FC with a large prior.
up.counts <- assay(out, "up")
down.counts <- assay(out, "down")
dir.logFC <- log2((up.counts+10)/(down.counts+10))
dir.logFC

# Calculating directionality index with Dixon's method.
dixon.stat <- sign(up.counts-down.counts)*2*((up.counts-down.counts)/(up.counts+down.counts))^2
dixon.stat

\dontshow{
unlink(fout, recursive=TRUE)
}
}

