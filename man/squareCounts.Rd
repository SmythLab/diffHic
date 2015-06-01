\name{squareCounts}
\alias{squareCounts}

\title{Load Hi-C interaction counts}

\description{Collate count combinations for interactions between pairs of bins across multiple Hi-C libraries.}

\usage{
squareCounts(files, param, width=50000, filter=1L)
}

\arguments{
	\item{files}{a character vector containing paths to the index files generated from each Hi-C library}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{width}{an integer scalar specifying the width of each bin in base pairs}
	\item{filter}{an integer scalar specifying the minimum count for each square}
}

\value{
A \code{DIList} object is returned containing the number of read pairs for each bin pair across all libraries.
The genomic coordinates for each bin are also recorded.
}

\details{
The genome is first split into non-overlapping adjacent bins of size \code{width}.
These bins are rounded to the nearest restriction site, forming the sides of potential squares.
The number of restriction fragments in each bin is stored as \code{nfrags} in the metadata of the output \code{region}.
Each square represents an interaction between the corresponding rounded bins on the genome. 
The number of read pairs between each pair of sides is counted for each library to obtain the count for the corresponding square.

% Can't use sliding mats because they're too difficult to interpret. Unlike ChIP-seq, there is no easy aggregation as there are no
% obvious peaks. Abundance aggregation will generally give meaningless results as you end up with massive clusters for local interactions.
% Without clustering, that means that each mat position must be reported separately. This is troublesome if they're overlapping as 
% we end up with problems with liberalness when interpreting the FDR.

Larger counts can be collected by increasing the value of \code{width}. 
This can improve detection power by increasing the evidence for significant differences.
However, this comes at the cost of spatial resolution as adjacent events in the same bin or square can no longer be distinguished.
This may reduce detection power if counts for differential interactions are contaminated by counts for non-differential interactions.

Low abundance squares with count sums below \code{filter} are not reported.
This reduces memory usage for large datasets.
These squares are probably uninteresting as detection power will be poor for low counts.
Another option is to increase \code{width} to reduce the total number of bins in the genome (and hence, the possible number of bin pairs).

% Self-circle PETs are artifacts resulting from inefficient cross-linking. As
% such, they are removed by default. The option to keep them is provided as some
% self-circles can form naturally e.g. in the genuine absence of interactions or with
% strong intra-fragment interactions. That said, they are difficult to interpret
% in the square-counting framework so removal is still recommended.
%
% It's hard to know which it is. Gets counted in the same squares, but can't tell if it's consistent with more or less compaction.

Counting will consider the values of \code{restrict}, \code{discard} and \code{cap} in \code{param}. 
See \code{\link{pairParam}} for more details.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, file=fout))

# Collating to count combinations.
y <- squareCounts(fout, param)
head(counts(y))
y <- squareCounts(fout, param, filter=1)
head(counts(y))
y <- squareCounts(fout, param, width=50, filter=1)
head(counts(y))
y <- squareCounts(fout, param, width=100, filter=1)
head(counts(y))

# Attempting with other parameters.
y <- squareCounts(fout, reform(param, restrict="chrA"), width=100, filter=1)
head(counts(y))
y <- squareCounts(fout, filter=1,
    param=reform(param, restrict=cbind("chrA", "chrB"))) 
head(counts(y))
y <- squareCounts(fout, filter=1,
    param=reform(param, cap=1), width=100)
head(counts(y))
y <- squareCounts(fout, width=100, filter=1, 
    param=reform(param, discard=GRanges("chrA", IRanges(1, 50))))
head(counts(y))

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\author{Aaron Lun}

\seealso{
\code{\link{preparePairs}},
\code{\link{cutGenome}},
\code{\link{DIList-class}}
}

\references{
Imakaev M et al. (2012). Iterative correction of Hi-C data reveals hallmarks of chromosome organization. \emph{Nat. Methods} 9, 999-1003.

Lieberman-Aiden E et al. (2009). Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome. \emph{Science} 326, 289-293. 
}

\keyword{counting}
