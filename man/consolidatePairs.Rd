\name{consolidatePairs}
\alias{consolidatePairs}

\title{Consolidate results for bin pairs}
\description{Consolidate differential testing results for bin pairs of differing sizes.}

\usage{
consolidatePairs(data.list, result.list, equiweight=TRUE, combine.args=list())
}

\arguments{
\item{data.list}{a list of \code{DIList} objects containing bin pairs}
\item{result.list}{a list of data frames containing the DB test results for each entry of \code{data.list}}
\item{equiweight}{a logical scalar indicating whether equal weighting from each bin size should be enforced}
\item{combine.args}{a list of parameters to pass to \code{\link{combineTests}}}
}

\value{
A list is returned, containing:
\item{id}{a list of integer vectors, where each vector corresponds to an object in \code{data.list}; the entries of the vector specify the parent bin pair to which each row of that object is assigned}
\item{pairs}{a \code{DIList} object containing the genomic coordinates of the parent bin pairs; see \code{\link{boxPairs}} for details}
\item{table}{a data frame containing the combined DB results for each \code{region}}
}

\details{
This function provides a convenient method to consolidate test results for bin pairs of differing sizes.
The interaction space is partitioned according to the largest bin size of all elements in \code{data.list}, to form parent bin pairs.
Each bin pair of each size in \code{data.list} is then assigned to the parent bin pair in which it is nested.
Note that the largest bin size must be an integer multiple of all smaller bin sizes; see \code{\link{boxPairs}} for more details.

For all nested bin pairs within a parent bin pair, the associated p-values are combined in \code{\link{combineTests}} using a weighted version of Simes' method.
This yields a single combined p-value, representing the evidence against the global null.
When \code{equiweight=TRUE}, the weight of a p-value of each bin pair is inversely proportional to the number of bin pairs of the same size in that parent bin pair.
This ensures that the results are not dominated by numerous smaller bin pairs.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{combineTests}},
\code{\link{boxPairs}}
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
out <- consolidatePairs(list(y, y2), list(result1, result2))
head(out$table)
consolidatePairs(list(y, y2), list(result1, result2), equiweight=FALSE)
head(out$table)
regions(out$pairs)

# Repeating with three sizes.
a10 <- a/10
b10 <- b/10
all.combos3 <- combn(length(cuts)/10, 2) # Bin size of 10.
y3 <- DIList(matrix(0, ncol(all.combos3), 1), anchors=all.combos3[2,],
	targets=all.combos3[1,], exptData=List(param=param, width=10),
	regions=GRanges(rep(c("chrA", "chrB"), c(a10, b10)), 
		IRanges(c((1:a10-1)*10+1, (1:b10-1)*10+1), c(1:a10*10, 1:b10*10))))
result3 <- data.frame(logFC=rnorm(nrow(y3)), PValue=runif(nrow(y3)), logCPM=0)

out <- consolidatePairs(list(y, y2, y3), list(result1, result2, result3))
head(out$table)
regions(out$pairs)
}

\keyword{testing}
