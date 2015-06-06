\name{consolidatePairs}
\alias{consolidatePairs}

\title{Consolidate results for interactions}
\description{Consolidate differential testing results for interactions from separate analyses.}

\usage{
consolidatePairs(indices, result.list, equiweight=TRUE, combine.args=list())
}

\arguments{
\item{indices}{a list of index vectors, specifying the cluster ID to which each interaction belongs}
\item{result.list}{a list of data frames containing the DB test results for each interaction}
\item{equiweight}{a logical scalar indicating whether equal weighting from each bin size should be enforced}
\item{combine.args}{a list of parameters to pass to \code{\link{combineTests}}}
}

\value{
A data frame is returned containing the combined DB results for each cluster.
}

\details{
Interactions from different analyses can be aggregated together using \code{\link{boxPairs}} or \code{\link{clusterPairs}}.
For example, test results can be consolidated for bin pairs of differing sizes.
This usually produces a \code{indices} vector that can be used as an input here.
Briefly, each vector in \code{indices} should correspond to one analysis, and each entry of that vector should correspond to an analyzed interaction.
The vector itself holds cluster IDs, such that interactions within/between analyses with the same ID belong in the same cluster.

For all bin pairs in a cluster, the associated p-values are combined in \code{\link{combineTests}} using a weighted version of Simes' method.
This yields a single combined p-value, representing the evidence against the global null.
When \code{equiweight=TRUE}, the weight of a p-value of each bin pair is inversely proportional to the number of bin pairs of the same size in that parent bin pair.
This ensures that the results are not dominated by numerous smaller bin pairs.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{combineTests}},
\code{\link{boxPairs}},
\code{\link{clusterPairs}}
}

\examples{
# Setting up the objects.
a <- 10
b <- 20
cuts <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))
param <- pairParam(cuts)

all.combos <- combn(length(cuts), 2) # Bin size of 1.
y <- DIList(matrix(0, ncol(all.combos), 1), anchors=all.combos[2,],
    targets=all.combos[1,], regions=cuts, exptData=List(param=param, width=1))

a5 <- a/5
b5 <- b/5
all.combos2 <- combn(length(cuts)/5, 2) # Bin size of 5.
y2 <- DIList(matrix(0, ncol(all.combos2), 1), anchors=all.combos2[2,],
    targets=all.combos2[1,], exptData=List(param=param, width=5),
	regions=GRanges(rep(c("chrA", "chrB"), c(a5, b5)), 
		IRanges(c((1:a5-1)*5+1, (1:b5-1)*5+1), c(1:a5*5, 1:b5*5))))

result1 <- data.frame(logFC=rnorm(nrow(y)), PValue=runif(nrow(y)), logCPM=0)
result2 <- data.frame(logFC=rnorm(nrow(y2)), PValue=runif(nrow(y2)), logCPM=0)

# Consolidating.
boxed <- boxPairs(y, y2)
out <- consolidatePairs(boxed$indices, list(result1, result2))
head(out)
out <- consolidatePairs(boxed$indices, list(result1, result2), equiweight=FALSE)
head(out)

# Repeating with three sizes.
a10 <- a/10
b10 <- b/10
all.combos3 <- combn(length(cuts)/10, 2) # Bin size of 10.
y3 <- DIList(matrix(0, ncol(all.combos3), 1), anchors=all.combos3[2,],
	targets=all.combos3[1,], exptData=List(param=param, width=10),
	regions=GRanges(rep(c("chrA", "chrB"), c(a10, b10)), 
		IRanges(c((1:a10-1)*10+1, (1:b10-1)*10+1), c(1:a10*10, 1:b10*10))))
result3 <- data.frame(logFC=rnorm(nrow(y3)), PValue=runif(nrow(y3)), logCPM=0)

boxed <- boxPairs(y, y2, y3)
out <- consolidatePairs(boxed$indices, list(result1, result2, result3))
head(out)
}

\keyword{testing}
