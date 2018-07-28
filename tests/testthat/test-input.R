# This tests non-standard methods for creating InteractionSet objects. 
# library(testthat); library(diffHic); source("test-input.R")

library(Matrix)
tmp.loc <- tempfile()
dir.create(tmp.loc)

###########################################################

# Mocking up some MTX and BED files.
set.seed(110000)
A <- rsparsematrix(1000, 1000, density=0.1, symmetric=TRUE, rand.x=function(n) round(runif(n, 1, 100)))
A.name <- file.path(tmp.loc, "A.mtx")
writeMM(file=A.name, A)

B <- rsparsematrix(1000, 1000, density=0.1, symmetric=TRUE, rand.x=function(n) round(runif(n, 1, 100)))
B.name <- file.path(tmp.loc, "B.mtx")
writeMM(file=B.name, B)

GR <- GRanges(sample(c("chrA", "chrB", "chrC"), 1000, replace=TRUE),
    IRanges(start=round(runif(1000, 1, 10000)),
        width=round(runif(1000, 50, 500))))
GR <- sort(GR)
bed.name <- file.path(tmp.loc, "regions.bed")
rtracklayer::export.bed(GR, con=bed.name)

test_that("readMTX2IntSet works as expected", {
    obs <- readMTX2IntSet(c(A.name, B.name), bed.name)
    expect_identical(regions(obs), GR)
    expect_type(assay(obs), "integer")
    expect_identical(assayNames(obs), "counts")

    # Checking the input against mergeCMs.
    cmA <- ContactMatrix(A, seq_len(1000), seq_len(1000), GR)
    cmB <- ContactMatrix(B, seq_len(1000), seq_len(1000), GR)
    ref <- mergeCMs(cmA, cmB)
    ref <- sort(ref)

    expect_identical(interactions(ref), interactions(obs))
    expect_equal(assay(ref), assay(obs))
    expect_identical(ref$totals, obs$totals)

    # Checking correct behaviour with only one file.        
    obs3 <- readMTX2IntSet(A.name, bed.name, as.integer=FALSE)
    ref <- mergeCMs(cmA)
    ref <- sort(ref)

    expect_identical(interactions(ref), interactions(obs3))
    expect_equal(assay(ref), assay(obs3))
    expect_identical(ref$totals, obs3$totals)

    # Checking for consistent behaviour when not dealing with integers.
    obs2 <- readMTX2IntSet(c(A.name, B.name), bed.name, as.integer=FALSE)
    expect_type(assay(obs2), "double")
    expect_equal(assay(obs2), assay(obs))
})

test_that("readMTX2IntSet behaves with silly inputs", {
    # No files.
    silly <- readMTX2IntSet(character(0), bed.name)
    expect_identical(dim(silly), c(0L, 0L))
    expect_identical(regions(silly), GR)

    # Empty files.
    A2 <- A
    A2[] <- 0
    A2.name <- file.path(tmp.loc, "A2.mtx")
    writeMM(file=A2.name, A2)
    silly <- readMTX2IntSet(A2.name, bed.name)
    expect_identical(dim(silly), c(0L, 1L))
    expect_identical(regions(silly), GR)

    # Inconsistent dimensions.
    A2 <- A[1:10,1:10]
    writeMM(file=A2.name, A2)
    expect_error(silly <- readMTX2IntSet(A2.name, bed.name), "not equal")
})

###########################################################

set.seed(110001)
chrs <- c(chrA=11, chrB=22, chrC=3)

test_that("mergeCMs works correctly", {
    for (nr in c(13, 27)) { 
    for (nc in c(13, 27)) { 
    for (lambda in c(1, 10)) { 
    for (mat.type in c("matrix", "dgCMatrix")) { 
    for (nsamples in 1:2) {
        # Simulating data under a range of scenarios.    
        N <- sum(chrs)
        all.starts <- round(runif(N, 1, 100))
        all.ends <- all.starts + round(runif(N, 5, 20))
        all.regions <- GRanges(rep(names(chrs), chrs), IRanges(all.starts, all.ends))
            
        all.anchor1 <- sample(N, nr)
        all.anchor2 <- sample(N, nc)
    
        collected <- list()
        for (i in seq_len(nsamples)) { 
            counts <- matrix(rpois(N*N, lambda=lambda), N, N)
            counts <- as.matrix(forceSymmetric(counts)) # Force symmetry, avoid considering non-identical redundant interactions.
            counts <- counts[all.anchor1, all.anchor2, drop=FALSE]
            counts <- as(counts, mat.type)
            x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)
            collected[[i]] <- x
        }

        # With no filtering.
        output <- do.call(mergeCMs, collected)
        for (i in seq_len(nsamples)) {
            current <- collected[[i]]
            test <- deflate(current, use.zero=FALSE)
            interactions(test) <- as(interactions(test), "ReverseStrictGInteractions")
            expect_identical(regions(test), regions(output))

            m <- match(test, output)
            expect_true(!any(is.na(m)))
            leftovers <- assay(output)[-m,i]
            expect_true(all(leftovers==0))

            ax <- anchors(test, id=TRUE)
            expect_identical(anchors(output, type="first", id=TRUE)[m], pmax(ax$first, ax$second))
            expect_identical(anchors(output, type="second", id=TRUE)[m], pmin(ax$first, ax$second))
            expect_identical(assay(test)[,1], assay(output)[m,i])
        }
    }}}}}
})

test_that("mergeCMs throws errors correctly", {
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
         
    expect_error(mergeCMs(x, x[,1]), "should be the same")
    expect_error(mergeCMs(x, x[1,]), "should be the same")
    x2 <- x
    regions(x2) <- resize(regions(x2), width=100)
    expect_error(mergeCMs(x, x2), "should be the same")
})
