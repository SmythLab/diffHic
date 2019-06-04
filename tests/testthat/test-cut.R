# This script is designed to test the restriction site finder. 
# library(testthat); library(diffHic); source("test-cut.R")

library(Biostrings)
library(BSgenome)

findRestrictionSites <- function(bs, pattern) {
	ps <- DNAString(pattern)
	remainder <- as.integer(nchar(pattern)/2 + 0.5001) # accounts for oddness.

	out <- list()
	for (chr in seqnames(bs)) {
		x <- matchPattern(pattern, bs[[chr]], fixed="subject")
		cuts <- start(x) + remainder - 1L
		out[[chr]] <- c(cuts, length(bs[[chr]]))
	}
	out
}

# We use the E.coli genome to do so 
# because it's short and we won't be spending ages putting it together.	
suppressPackageStartupMessages(require(BSgenome.Ecoli.NCBI.20080805))

test_that("single cutter applications of cutGenome work correctly", {
    comp <- function(target, overhang=4L) {
        current <- cutGenome(Ecoli, target, overhang=overhang)
        theoretical <- findRestrictionSites(Ecoli, target)
        expect_identical(sort(names(theoretical)), sort(seqlevels(current)))

        half.ohang <- as.integer(overhang/2L)
        for (x in names(theoretical)) {
            n <- length(theoretical[[x]])
            cuts1 <- end(current[x==seqnames(current)])
            cuts1[-n] <- cuts1[-n]-half.ohang
            expect_identical(cuts1, theoretical[[x]])
        }
        return(head(current))
    }

    comp("GGATCC", 4L) # BamHI

    comp("GAATTC", 4L) # EcoRI

    comp("AACGTT", 2L) # AclI

    comp("GCGC", 2L) # HinP1I

    comp("GATC", 4L) # MboI

    comp("GCGGCCGC", 4L) # NotI

    comp("GGCGCGCC", 4L) # AscI

    comp("ACNGT", 3L) # infI
})

test_that("double cutter applications of cutGenome work correctly", {
    bamhI <- cutGenome(Ecoli, "GGATCC", 4L)
    aclI <- cutGenome(Ecoli, "AACGTT", 2L)
    both <- cutGenome(Ecoli, c("GGATCC", "AACGTT"), c(4L, 2L))

    # Manually creating 'both'.
    olap <- findOverlaps(aclI, bamhI)
    combined <- pintersect(aclI[queryHits(olap)], bamhI[subjectHits(olap)])
    combined$hit <- NULL

    expect_identical(both, combined)
    expect_identical(sort(combined), combined)
})

test_that("cutGenome throws errors correctly", {
    expect_error(cutGenome(Ecoli, "GATC", 1L), "even")
    expect_error(cutGenome(Ecoli, "GATC", 6L), "non-negative")
    expect_error(cutGenome(Ecoli, "GATC", -1L), "non-negative")
    expect_error(cutGenome(Ecoli, "GCTC", -1L), "inverse palindrome")
})

