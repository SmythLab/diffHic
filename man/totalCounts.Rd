\name{totalCounts}
\alias{totalCounts}

\title{Get the total counts}

\description{Get the total number of read pairs in a set of Hi-C libraries.}

\usage{
totalCounts(files, param)
}

\arguments{
	\item{files}{a character vector containing paths to the index files generated from each Hi-C library}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
}

\value{
An integer vector is returned containing the total number of read pairs in each library.
}

\details{
As the name suggests, this function counts the total number of read pairs in each index file prepared by \code{\link{preparePairs}}.
Use of \code{param$fragments} ensures that the chromosome names in each index file are consistent with those in the desired genome (e.g., from \code{\link{cutGenome}}).
Counting will also consider the values of \code{restrict}, \code{discard} and \code{cap} in \code{param}.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, file=fout))

# Counting totals, and comparing them.
totalCounts(fout, param)
squareCounts(fout, param, width=10)$totals

new.param <- reform(param, restrict="chrA")
totalCounts(fout, new.param)
squareCounts(fout, new.param, width=10)$totals

new.param <- reform(param, discard=GRanges("chrA", IRanges(1, 50)))
totalCounts(fout, new.param)
squareCounts(fout, new.param, width=10)$totals

new.param <- reform(param, cap=1)
totalCounts(fout, new.param)
squareCounts(fout, new.param, width=10)$totals

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\author{Aaron Lun}

\seealso{
\code{\link{preparePairs}}, 
\code{\link{cutGenome}},
\code{\link{pairParam}},
\code{\link{squareCounts}}
}

\keyword{counting}
