\name{plotDI}
\alias{plotDI}
\alias{rotDI}

\title{Construct a plaid plot of differential interactions}

\description{Plot differential interactions in a plaid format with informative colouring.}

\usage{
plotDI(data, fc, anchor, target=anchor, col.up="red", col.down="blue",
    background="grey70", zlim=NULL, xlab=NULL, ylab=NULL, diag=TRUE, ...)
rotDI(data, fc, region, col.up="red", col.down="blue",
    background="grey70", zlim=NULL, xlab=NULL, ylab="Gap", ...)

}

\arguments{
	\item{data}{a \code{DIList} object}
	\item{fc}{a numeric vector of log-fold changes}
	\item{anchor}{a \code{GRanges} object of length 1 specifying the anchor region}
	\item{target}{a \code{GRanges} object of length 1 specifying the target region}
	\item{region}{a \code{GRanges} object of length 1 specifying the region of interest}
	\item{col.up}{any type of R color to describe the maximum color for positive log-fold changes}
	\item{col.down}{any type of R color to describe the maximum color for negative log-fold changes}
	\item{background}{any type of R color, specifying the background color of the interaction space}
	\item{zlim}{a numeric scalar indicating the maximum absolute log-fold change}
	\item{xlab}{character string for the x-axis label on the plot, defaults to the anchor chromosome name}
	\item{ylab}{character string for the y-axis label on the plot, defaults to the target chromosome name in \code{plotDI}}
	\item{diag}{a logical scalar specifying whether boxes should be shown above the diagonal for intra-chromosomal plots in \code{plotDI}}
	\item{...}{other named arguments to be passed to \code{\link{plot}}}
}	

\details{
The \code{plotDI} function constructs a plaid plot on the current graphics device.
The \code{anchor} and \code{target} regions are represented by the x- and y-axes, respectively.
Each bin pair is represented by a box in the plotting space, where each side of the box represents a bin.
Plotting space that is not covered by any bin pair is shown in \code{background}.

The colour of the box depends on the magnitude and sign of the log-fold change in \code{fc}.
Positive log-FCs will range from white to \code{col.up}, whereas negative log-FCs will range from white to \code{col.down}.
The chosen color is proportional to the magnitude of the log-FC, and the most extreme colors are only obtained at the maximum absolute log-FC in \code{fc}.
The maximum value can be capped at \code{zlim} for better resolution of small log-FCs.

If \code{diag=TRUE}, boxes will also be plotted above the diagonal for intra-chromosomal plots. 
This is set as the default to avoid confusion when anchor/target ordering is misspecified, but can be turned off to reduce redundancy in visualization around the diagonal.

The \code{rotDI} function constructs a rotated plot of differential interactions, for visualization of local changes.
See \code{\link{rotPlaid}} for more details. 
}

\value{
A (rotated) plaid plot of differential interactions is produced on the current graphics device. 
A function is also invisibly returned that converts log-FCs into colors.
This is useful for coordinating the colors, e.g., when constructing a separate color bar.
}

\author{Aaron Lun}

\examples{
# Setting up the objects.
a <- 10
b <- 20
regions <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)),
    seqinfo=Seqinfo(seqlengths=c(chrA=a, chrB=b), seqnames=c("chrA", "chrB")))

set.seed(3423)
all.anchors <- sample(length(regions), 500, replace=TRUE)
all.targets <- as.integer(runif(500, 1, all.anchors+1))
out <- DIList(matrix(0, 500, 1), anchors=all.anchors, targets=all.targets, 
    regions=regions, exptData=List(width=1))
fc <- runif(nrow(out), -2, 2)

# Constructing intra-chromosomal DI plots around various regions
plotDI(out, fc, anchor=GRanges("chrA", IRanges(1, 10)), 
    target=GRanges("chrA", IRanges(1, 10)), diag=TRUE)
plotDI(out, fc, anchor=GRanges("chrA", IRanges(1, 10)), 
    target=GRanges("chrA", IRanges(1, 10)), diag=FALSE)

# Constructing inter-chromosomal DI plots around various regions
xxx <- plotDI(out, fc, anchor=GRanges("chrB", IRanges(1, 10)), 
    target=GRanges("chrA", IRanges(1, 20)), diag=TRUE)
plotDI(out, fc, anchor=GRanges("chrB", IRanges(1, 100)), 
    target=GRanges("chrA", IRanges(1, 200)), diag=TRUE, zlim=5)

# Making colorbars.
xxx((-10):10/10)
xxx((-20):20/20)

# Rotated.
rotDI(out, fc, region=GRanges("chrA", IRanges(1, 200)))
rotDI(out, fc, region=GRanges("chrB", IRanges(1, 200)))
}

\seealso{
	\code{\link{plotPlaid}},
	\code{\link{rotPlaid}},
	\code{\link{squareCounts}}
}

\references{
Lieberman-Aiden E et al. (2009). Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome. \emph{Science} 326, 289-293. 
}

\keyword{visualization}
