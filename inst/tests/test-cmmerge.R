####################################################################################################
# Checking the mergeCMs function.

suppressPackageStartupMessages(library(diffHic))
suppressPackageStartupMessages(library(Matrix))

cmcomp <- function(chrs, Nr, Nc, nsamples, lambda=10, mat.type="matrix", filter=10) {
    N <- sum(chrs)
    all.starts <- round(runif(N, 1, 100))
    all.ends <- all.starts + round(runif(N, 5, 20))
    all.regions <- GRanges(rep(names(chrs), chrs), IRanges(all.starts, all.ends))
        
    all.anchor1 <- sample(N, Nr)
    all.anchor2 <- sample(N, Nc)

    original <- collected <- list()
    for (i in seq_len(nsamples)) { 
        counts <- matrix(rpois(N*N, lambda=lambda), N, N)
        counts <- as.matrix(forceSymmetric(counts)) # Force symmetry, avoid considering non-identical redundant interactions.
        counts <- counts[all.anchor1, all.anchor2, drop=FALSE]
        original[[i]] <- counts

        counts <- as(counts, mat.type)
        x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)
        collected[[i]] <- x
    }

    # Checking that this is the same as sample-wise inflations.
    output <- do.call(mergeCMs, collected)
    for (i in seq_along(nsamples)) {
        current <- collected[[i]]
        test <- deflate(current, use.zero=TRUE) # Forcing use of zeros for sparse matrices.
        stopifnot(identical(regions(test), regions(output)))

        ax <- anchors(test, id=TRUE)
        stopifnot(identical(anchors(output, type="first", id=TRUE), pmax(ax$first, ax$second)))
        stopifnot(identical(anchors(output, type="second", id=TRUE), pmin(ax$first, ax$second)))
        stopifnot(identical(assay(test)[,1], assay(output)[,i]))
    }

    # Checking that the totals are the same.
    row.space <- as.integer(matrix(all.anchor1, Nr, Nc))
    col.space <- as.integer(matrix(all.anchor2, Nr, Nc, byrow=TRUE))
    combined <- cbind(pmax(col.space, row.space), pmin(col.space, row.space))
    keep <- !duplicated(combined)
    stopifnot(isTRUE(all.equal(output$totals, sapply(original, FUN=function(M) { sum(M[keep]) })))) 

    return(head(assay(output)))
}

set.seed(10000)
chrs <- c(chrA=10, chrB=10)
cmcomp(chrs, 10, 10, nsamples=2, lambda=5)
cmcomp(chrs, 10, 10, nsamples=4, lambda=2)
cmcomp(chrs, 15, 5, nsamples=2, lambda=5)
cmcomp(chrs, 15, 5, nsamples=4, lambda=2)
cmcomp(chrs, 5, 15, nsamples=2, lambda=5)
cmcomp(chrs, 5, 15, nsamples=4, lambda=2)

# Testing out behaviour with sparse matrices.
cmcomp(chrs, 10, 10, nsamples=2, lambda=5, mat.type="dgCMatrix")
cmcomp(chrs, 10, 10, nsamples=4, lambda=2, mat.type="dgCMatrix")
cmcomp(chrs, 10, 10, nsamples=10, lambda=1, mat.type="dgCMatrix")
cmcomp(chrs, 10, 10, nsamples=2, lambda=5, mat.type="dgeMatrix")
cmcomp(chrs, 10, 10, nsamples=4, lambda=2, mat.type="dgeMatrix")

# Testing out behaviour with more chromosomes.
chrs <- c(chrA=15, chrB=10, chrC=5)
cmcomp(chrs, 10, 10, nsamples=2, lambda=5)
cmcomp(chrs, 10, 10, nsamples=4, lambda=2)
cmcomp(chrs, 15, 5, nsamples=2, lambda=5)
cmcomp(chrs, 15, 5, nsamples=4, lambda=2)
cmcomp(chrs, 5, 15, nsamples=2, lambda=5)
cmcomp(chrs, 5, 15, nsamples=4, lambda=2)

# Ensuring that an error is chucked when the regions involved are not identical.
    
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Nr <- 10
Nc <- 20
all.anchor1 <- sample(N, Nr)
all.anchor2 <- sample(N, Nc)
counts <- matrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)
     
try(mergeCMs(x, x[,1]))
try(mergeCMs(x, x[1,]))
x2 <- x
regions(x2) <- resize(regions(x2), width=100)
try(mergeCMs(x, x2))

##################################################################################################

