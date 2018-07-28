# This tests non-standard methods for creating InteractionSet objects. 
# library(testthat); library(diffHic); source("test-input.R")

library(Matrix)
tmp.loc <- tempfile()
dir.create(tmp.loc)

###########################################################

# Mocking up some MTX and BED files.
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
    # Checking the input against mergeCMs.
    obs <- readMTX2IntSet(c(A.name, B.name), bed.name)
    expect_identical(regions(obs), GR)

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
})

test_that("readMTX2IntSet options work as expected", {
    obs <- readMTX2IntSet(c(A.name, B.name), bed.name)
    expect_type(assay(obs), "integer")
    expect_identical(assayNames(obs), "counts")

    # Checking for consistent behaviour when not dealing with integers.
    obs2 <- readMTX2IntSet(c(A.name, B.name), bed.name, as.integer=FALSE)
    expect_type(assay(obs2), "double")
    expect_equal(assay(obs2), assay(obs))

    # Responds to alternative names for consistent behaviour when not dealing with integers.
    obs2 <- readMTX2IntSet(c(A.name, B.name), bed.name, assay.type="whee")
    expect_identical(assayNames(obs2), "whee")
    expect_identical(assay(obs), assay(obs2))
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
