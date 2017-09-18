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
head(pairdist(data))

anchors(data, type="first")
anchors(data, type="second")
assay(data)
regions(data)

data$totals
colData(data)
metadata(data)

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
dummy <- InteractionSet(matrix(as.integer(rpois(npairs*nlibs, runif(npairs, 10, 100))), nrow=npairs),
    colData=DataFrame(totals=runif(nlibs, 1e6, 2e6)), 
    GInteractions(anchor1=anchors, anchor2=targets, regions=GRanges("chrA", IRanges(1:npts, 1:npts)), mode="reverse"))

normOffsets(dummy, se.out=FALSE)
normOffsets(dummy, logratioTrim=0, se.out=FALSE)
normOffsets(dummy, sumTrim=0.2, se.out=FALSE)
head(normOffsets(dummy, type="loess", se.out=FALSE))
head(normOffsets(dummy, type="loess", span=0.5, se.out=FALSE))

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
