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
For each bin, this function computes the total count for all bin pairs involving that bin and the \code{span} upstream bins.
This is repeated for the \code{span} downstream bins.
If multiple \code{files} are specified, the total of the averages across all libraries for each bin pair is computed instead.

The total up- and downstream counts are returned in the \code{Up} and \code{Down} fields of the output \code{GRanges} object.
This can be used to compute a directionality statistic, e.g., as defined by Dixon et al, or by computing the log-fold change between fields.
The latter is easier to compute and model with a Gaussian distribution.
In any case, the directionalities can be used in a HMM to identify domains - see the user's guide for more details.
}

\value{
A \code{GRanges} object with one entry for each bin in the genome.
The \code{Up} and \code{Down} fields in the metadata contain the total average count for the up-/downstream interactions involving each bin.
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
domainDirections(fout, param, width=10, span=1)

# Calculating directionality log-FC with a large prior.
dir.logFC <- log2((out$Up+10)/(out$Down+10))
dir.logFC

# Calculating directionality index with Dixon's method.
dixon.stat <- sign(out$Up-out$Down)*2*((out$Up-out$Down)/(out$Up+out$Down))^2
dixon.stat

\dontshow{
unlink(fout, recursive=TRUE)
}
}

