#################################################################
# This tests the filterDirect machinery.

require(diffHic)
set.seed(1212)

ends <- c(1:10, 1:5)*10L
starts <- ends - 9L
my.regions <- GRanges(rep(c("chrA", "chrB"), c(10, 5)), IRanges(starts, ends))

all.inters <- expand.grid(1:10, 11:15) # i.e., only inter-chromosomal interactions here.
x <- InteractionSet(matrix(10, nrow=nrow(all.inters), 1), 
                    GInteractions(all.inters[,1], all.inters[,2], my.regions, mode="reverse"),
                    colData=DataFrame(totals=1e6)) 

# Vanilla comparison.

out <- filterDirect(x)
stopifnot(length(out$abundances) && length(out$threshold))
stopifnot(all(out$abundances==out$threshold))
out$threshold

out <- filterDirect(x, prior.count=5)
stopifnot(length(out$abundances) && length(out$threshold))
stopifnot(all(out$abundances==out$threshold))
out$threshold

# Comparison with some missing entries; but these should not be the median.

out <- filterDirect(x)

x2 <- x
assay(x2)[1:20,] <- 50
out2 <- filterDirect(x2)
stopifnot(out2$threshold==out$threshold)
stopifnot(out2$threshold==median(out2$abundances))
out2$threshold

x2 <- x2[1:30,]
out2 <- filterDirect(x2)
stopifnot(out2$threshold==out$threshold)
stopifnot(out2$threshold!=median(out2$abundances))
out2$threshold

x2 <- x2[1,]
out2 <- filterDirect(x2) # Checking makeEmpty()
stopifnot(out2$threshold==edgeR::aveLogCPM(0, 1e6))

# Seeing what happens if we add in a whole bunch of intra-chromosomals.

all.combos <- combn(10, 2)
xi <- InteractionSet(matrix(100, nrow=ncol(all.combos), 1), 
                    GInteractions(all.combos[1,], all.combos[2,], my.regions, mode="reverse"),
                    colData=DataFrame(totals=1e6)) 
xi <- rbind(x, xi, xi)
outi <- filterDirect(xi)
stopifnot(outi$threshold==out$threshold)
stopifnot(outi$threshold!=median(outi$abundances))
outi$threshold

# Comparison with a reference object.

ends <- c(1:10, 1:5)*20L
starts <- ends - 19L
ref.regions <- GRanges(rep(c("chrA", "chrB"), c(10, 5)), IRanges(starts, ends))

all.inters <- expand.grid(1:10, 11:15)
ref <- InteractionSet(matrix(40, nrow=nrow(all.inters), 1), 
                    GInteractions(all.inters[,1], all.inters[,2], ref.regions, mode="reverse"),
                    colData=DataFrame(totals=1e6)) 

outr <- filterDirect(x, ref=ref)
stopifnot(all(outr$threshold==outr$abundances))
stopifnot(all(outr$threshold==outr$ref$abundances))
outr$threshold

outr <- filterDirect(x, ref=ref, prior.count=5)
stopifnot(all(outr$threshold==outr$abundances))
stopifnot(all(outr$threshold==outr$ref$abundances))
outr$threshold

subref <- ref
assay(subref)[1:20,] <- 100
out.sub <- filterDirect(x, ref=subref)
stopifnot(out.sub$threshold==out$threshold)
stopifnot(all(out.sub$abundances==out$abundances))
stopifnot(all(out.sub$ref$threshold==median(out.sub$ref$abundances)))
out.sub$threshold

subref <- subref[1:30,]
out.sub <- filterDirect(x, ref=subref)
stopifnot(out.sub$threshold==out$threshold)
stopifnot(all(out.sub$abundances==out$abundances))
stopifnot(all(out.sub$ref$threshold!=median(out.sub$ref$abundances)))
out.sub$threshold

refi <- InteractionSet(matrix(100, nrow=ncol(all.combos), 1), 
                    GInteractions(all.combos[1,], all.combos[2,], ref.regions, mode="reverse"),
                    colData=DataFrame(totals=1e6)) 
refi <- rbind(ref, refi, refi)
outi <- filterDirect(xi, ref=refi)
stopifnot(outi$threshold==out$threshold)
stopifnot(outi$threshold!=median(outi$abundances))
stopifnot(outi$threshold!=median(outi$ref$abundances))
outi$threshold

#################################################################
# This tests the filterTrended machinery, which is substantially more complex.

ends <- c(1:10, 1:5)*10L
starts <- ends - 9L
my.regions <- GRanges(rep(c("chrA", "chrB"), c(10, 5)), IRanges(starts, ends))

all.A1 <- sample(length(my.regions), 100, replace=TRUE)
all.A2 <- sample(length(my.regions), 100, replace=TRUE)
x <- InteractionSet(matrix(100, nrow=100, 1), 
                    GInteractions(all.A1, all.A2, my.regions, mode="reverse"),
                    colData=DataFrame(totals=1e6)) 

out <- filterTrended(x)
stopifnot(all.equal(out$log.distance, log10(pairdist(x)+median(width(my.regions)))))
head(out$log.distance)
stopifnot(all(out$abundances==edgeR::aveLogCPM(100, 1e6)))
head(out$abundances)

# Checking that we can fill in the missing distances.

is.intra <- !is.na(pairdist(x))
a.pts <- anchors(x, type="first", id=TRUE)[is.intra]
t.pts <- anchors(x, type="second", id=TRUE)[is.intra]

o <- order(a.pts, t.pts)
a.pts <- a.pts[o]
t.pts <- t.pts[o]

all.chrs <- seqnames(regions(x))
all.mids <- (start(regions(x))+end(regions(x)))/2
extra.dist <- .Call(diffHic:::cxx_get_missing_dist, cumsum(runLength(all.chrs)),
                    a.pts-1L, t.pts-1L, all.mids)

suppressWarnings(cm <- inflate(x, rows=NULL, columns=NULL)) # Reference way to do it.
ref.dist <- pairdist(cm)
ref.dist <- ref.dist[as.matrix(is.na(as.matrix(cm))) & !is.na(ref.dist) & lower.tri(ref.dist, diag=TRUE)]
stopifnot(all.equal(sort(extra.dist), sort(ref.dist)))
head(extra.dist)

# Checking that we fit a sensible trend.

fit <- limma::loessFit(c(out$abundances, rep(edgeR::aveLogCPM(0, 1e6), length(ref.dist))),
                       x=c(out$log.distance, log10(ref.dist + median(width(my.regions)))),
                       span=formals(filterTrended)$span)$fitted

ref.threshold <- fit[seq_along(out$abundances)]
dout <- filterDirect(x)
ref.threshold[is.na(ref.threshold)] <- dout$threshold
stopifnot(all.equal(ref.threshold, out$threshold))
head(out$threshold)

# Checking that behaviour upon setting 'ref' is reasonable.

xr <- x
assay(xr)[] <- 40
regions(xr) <- resize(regions(xr), width=20, fix="center")

outr <- filterTrended(x, ref=xr)
new.threshold <- approx(outr$ref$log.distance, outr$ref$threshold, xout=outr$log.distance, rule=2)$y
doutr <- filterDirect(x, ref=xr)
doutr$threshold

new.threshold[is.na(outr$log.distance)] <- doutr$threshold
stopifnot(all.equal(new.threshold, outr$threshold))
head(outr$threshold)

#################################################################
# End.

