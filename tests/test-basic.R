# These are just placeholders for the real things in inst/tests.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(fragments=cuts)

# Setting up the parameters
fout <- "output.h5"
preparePairs(hic.file, param, file=fout)
head(getPairData(fout, param))

loadChromos(fout)
head(loadData(fout, "chrA", "chrA"))
head(loadData(fout, "chrA", "chrB"))

# Loading the counts.
data <- squareCounts(fout, param, width=50, filter=1)
data

margins <- marginCounts(fout, param, width=50)
margins
totalCounts(fout, param)

regions <- GRanges("chrA", IRanges(c(1, 100, 150), c(20, 140, 160)))
connectCounts(fout, param, regions=regions, filter=1L)

# Checking some values.
head(getArea(data))
head(getDistance(data))

anchors(data)
targets(data)
counts(data)
regions(data)

data$totals
colData(data)
exptData(data)

asDGEList(data)
asDGEList(data, lib.size=20)$samples
asDGEList(data, norm.factors=2, group="a")$samples

# Simple normalization with dummy data.
set.seed(3423746)
npts <- 100
npairs <- 5000
nlibs <- 4
anchors <- sample(npts, npairs, replace=TRUE)
targets <- sample(npts, npairs, replace=TRUE)
dummy <- DIList(counts=matrix(rpois(npairs*nlibs, runif(npairs, 10, 100)), nrow=npairs),
    totals=runif(nlibs, 1e6, 2e6), anchors=pmax(anchors, targets), targets=pmin(anchors, targets),
    regions=GRanges("chrA", IRanges(1:npts, 1:npts)))

normalize(dummy)
normalize(dummy, logratio=0)
normalize(dummy, lib.sizes=c(10, 20, 15, 25))
head(normalize(dummy, type="loess"))
head(normalize(dummy, type="loess", span=0.5))

# Playing around with some bin counts.
stuff <- correctedContact(data)
head(stuff$truth)

data.large <- squareCounts(fout, param, width=100, filter=1)
boxed <- boxPairs(larger=data.large, smaller=data)
head(boxed$indices$larger)
head(boxed$indices$smaller)

head(enrichedPairs(data))
head(clusterPairs(data, tol=10)$indices[[1]])

# End.

unlink(fout)
