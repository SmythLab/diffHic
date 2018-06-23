\name{diClusters}
\alias{diClusters}

\title{Cluster significant bin pairs to DIs}
\description{Cluster significant bin pairs to DIs with post-hoc cluster-level FDR control.}

\usage{
diClusters(data.list, result.list, target, equiweight=TRUE, cluster.args=list(), 
    pval.col="PValue", fc.col=NA, grid.length=21, iterations=4)
}

\arguments{
\item{data.list}{an InteractionSet or a list of InteractionSet objects containing bin pairs.}
\item{result.list}{a data frame or a list of data frames containing the DI test results for each bin pair.}
\item{target}{a numeric scalar specifying the desired cluster-level FDR.}
\item{equiweight}{a logical scalar indicating whether equal weighting from each input object should be enforced}
\item{cluster.args}{a list of parameters to supply to \code{\link{clusterPairs}}.}
\item{pval.col}{a character string or integer scalar specifying the column of p-values for elements in \code{result.list}.}
\item{fc.col}{a character string or integer scalar specifying the column of log-fold changes for elements in \code{result.list}.}
\item{grid.length, iterations}{Parameters to supply to \code{\link{controlClusterFDR}}.}
}

\value{
A list of cluster indices and minimum bounding boxes is returned as described in \code{\link{clusterPairs}}.
An additional \code{FDR} field is also present, containing the estimate of the cluster-level FDR.
}

\details{
Bin pairs are identified as being significant based on the adjusted p-values in the corresponding data frame of \code{result.list}.
Only these significant bin pairs are clustered together via \code{\link{clusterPairs}}.
This identifies DIs consisting only of significant bin pairs.
By default, the \code{tol} parameter in \code{\link{clusterPairs}} is set to 1 bp, i.e., all adjacent bin pairs are clustered together.
If \code{fc.col} is specified, all clusters consist of bin pairs that are changing in the same direction.

The aim is to avoid very large clusters from blind clustering in areas of the interaction space that have high interaction intensity.
This includes interactions within structural domains, or in data sets where interactions are difficult to define due to high levels of noise.
Post-hoc control of the cluster-level FDR is performed using the \code{\link{controlClusterFDR}} function.
This is necessary as clustering is not blind to the test results.
By default, the cluster-level FDR is controlled at 0.05 if \code{target} is not specified.

Some effort is required to equalize the contribution of the results from each element of \code{result.list}.
This is done by setting \code{equiweight=TRUE}, where the weight of each bin pair is inversely proportional to the number of bin pairs from that analysis.
These weights are used as frequency weights for bin pair-level FDR control, when identifying significant bin pairs prior to clustering.
Otherwise, the final results would be dominated by large number of small bin pairs.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{clusterPairs}},
\code{\link{controlClusterFDR}}
}

\examples{
# Setting up the objects.
a <- 10
b <- 20
cuts <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))
param <- pairParam(cuts)

all.combos <- combn(length(cuts), 2) # Bin size of 1.
y1 <- InteractionSet(matrix(0, ncol(all.combos), 1), 
    GInteractions(anchor1=all.combos[2,], anchor2=all.combos[1,], regions=cuts, mode="reverse"),
    colData=DataFrame(lib.size=1000), metadata=List(param=param, width=1))

set.seed(1000)
result1 <- data.frame(logFC=rnorm(nrow(y1)), PValue=runif(nrow(y1)), logCPM=0)
result1$PValue[sample(nrow(result1), 50)] <- 0

# Consolidating with post-hoc control.
out <- diClusters(y1, result1, target=0.05, cluster.args=list(tol=1))
out
}

\keyword{testing}
