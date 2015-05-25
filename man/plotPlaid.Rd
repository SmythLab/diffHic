\name{plotPlaid}
\alias{plotPlaid}
\alias{rotPlaid}

\title{Construct a plaid plot of interactions}

\description{Plot interactions between two sequences in a plaid format with informative colouring.}

\usage{
plotPlaid(file, param, first.region, second.region=first.region, 
    width=10000, col="black", max.count=20, xlab=NULL, ylab=NULL, 
	diag=TRUE, count=FALSE, count.args=list(), ...)
rotPlaid(file, param, region, width=10000, col="black", 
    max.count=20, xlab=NULL, ylab="Gap", ...) 
}

\arguments{
	\item{file}{character string specifying the path to an index file produced by \code{\link{preparePairs}}}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{first.region}{a \code{GRanges} object of length 1 specifying the first region}
	\item{second.region}{a \code{GRanges} object of length 1 specifying the second region}
	\item{region}{a \code{GRanges} object of length 1 specifying the region of interest}
	\item{width}{an integer scalar specifying the width of each bin in base pairs}
	\item{col}{any type of R color to describe the color of the plot elements}
	\item{max.count}{a numeric scalar specifying the count for which the darkest color is obtained}
	\item{xlab}{character string for the x-axis label on the plot, defaults to the chromosome name of \code{first.region}}
	\item{ylab}{character string for the y-axis label on the plot, defaults to the chromosome name of \code{second.regeion} in \code{plotPlaid}}
	\item{diag}{a logical scalar specifying whether boxes should be shown above the diagonal for intra-chromosomal plots in \code{plotPlaid}}
	\item{count}{a logical scalar specifying whether the count for each bin should be plotted in \code{plotPlaid}}
	\item{count.args}{a named list of arguments to be passed to \code{\link{text}} for plotting of bin counts, if \code{count=TRUE}}
	\item{...}{other named arguments to be passed to \code{\link{plot}}}
}	

\details{
The \code{plotPlaid} function constructs a plaid plot on the current graphics device.
The intervals of the \code{first.region} and \code{second.region} are represented by the x- and y-axes, respectively.
Each region is partitioned into bins of size \code{width}.
Each bin pair is represented by a box in the plotting space, where each side of the box represents a bin.
The colour of the box depends on the number of read pairs mapped between the corresponding bins.

The resolution of colors can be controlled by varying \code{max.count}. 
All boxes with counts above \code{max.count} will be assigned the maximum intensity.
Other boxes will be assigned a color of intensity proportional to the size of the count, such that a count of zero results in white space.
Smaller values of \code{max.count} will improve contrast at low counts at the cost of contrast at higher counts. 
Scaling \code{max.count} is recommended for valid comparisons between libraries of different sizes (e.g., larger \code{max.count} for larger libraries).

If \code{count=TRUE}, the number of read pairs will be shown on top of each bin. 
This will be slower to plot but can be useful in some cases, e.g., when more detail is required, or when the range of colors is not sufficient to capture the range of counts in the data. 
If \code{diag=TRUE}, boxes will also be plotted above the diagonal for intra-chromosomal plots. 
This is set as the default to avoid confusion when \code{first.region} is not set as the anchor range, i.e., it has a lower sorting order than \code{second.region}.
However, this can also be turned off to reduce redundancy in visualization around the diagonal.

The \code{rotPlaid} function constructs a plaid plot that has been rotated by 45 degrees.
This is useful for visualizing local interactions within a specified \code{region}.
In a rotated plot, the x-coordinate of a box in the plotting space represents the midpoint between two interacting bins, while the y-coordinate represents the distance between bins.
More simply, the interacting bins of a box can be identified by tracing diagonals from the edges of the box to the x-axis.

Note that the plotted boxes for the bin pairs may overwrite the bounding box of the plot.
This can be fixed by running \code{box()} after each \code{plotPlaid} call.
}

\value{
A (rotated) plaid plot is produced on the current graphics device. 
For both functions, a function is invisibly returned that converts counts into colors.
This is useful for coordinating the colors, e.g., when constructing a separate color bar.
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
originals <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(originals)

# Setting up parameters
fout <- "temp_saved.h5"
invisible(preparePairs(hic.file, param, fout))

# Constructing intra-chromosomal plaid plots around various regions.
plotPlaid(fout, param, first.region=GRanges("chrA", IRanges(1, 100)), 
    second.region=GRanges("chrA", IRanges(1, 200)), width=50, diag=TRUE)
box()
xxx <- plotPlaid(fout, param, first.region=GRanges("chrA", IRanges(1, 100)), 
    second.region=GRanges("chrA", IRanges(1, 200)), width=50, diag=FALSE)

# Making colorbars.
xxx(1:2)
xxx(1:5)
xxx(1:10)

# Constructing inter-chromosomal plaid plots around various regions
plotPlaid(fout, param, first.region=GRanges("chrB", IRanges(1, 100)), 
    second.region=GRanges("chrA", IRanges(1, 200)), width=50)
plotPlaid(fout, param, first.region=GRanges("chrB", IRanges(1, 100)), 
    second.region=GRanges("chrA", IRanges(1, 200)), width=100)

# For a hypothetical second library which is half the size of the previous one:
plotPlaid(fout, param, first.region=GRanges("chrB", IRanges(1, 100)),
    second.region=GRanges("chrA", IRanges(1, 200)), width=100, max.count=20, count=TRUE)
plotPlaid(fout, param, first.region=GRanges("chrB", IRanges(1, 100)),
    second.region=GRanges("chrA", IRanges(1, 200)), width=100, max.count=40,
    count=TRUE, count.args=list(col="blue"))

# Rotated
rotPlaid(fout, param, region=GRanges("chrA", IRanges(1, 200)), width=50)
rotPlaid(fout, param, region=GRanges("chrA", IRanges(1, 200)), width=100)

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\seealso{
	\code{\link{preparePairs}}
}

\references{
Lieberman-Aiden E et al. (2009). Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome. \emph{Science} 326, 289-293. 
}

\keyword{visualization}
