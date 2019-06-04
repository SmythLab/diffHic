#' @export
#' @importFrom Biostrings DNAString reverseComplement readDNAStringSet matchPattern
#' @importClassesFrom BSgenome BSgenome
#' @importFrom methods is
#' @importFrom stats start
#' @importFrom GenomeInfoDb seqnames genome seqlevels<- seqinfo<- seqlengths Seqinfo
#' @importFrom BiocGenerics sort
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlaps pintersect
cutGenome <- function(bs, pattern, overhang=4L) 
# This finds the target cut sites in the genome. It currently only searches the
# sense strand, which is fine because if the patterns is an inverse palindrome.
# Otherwise, there may be some problems as you'd need to search the reverse
# strand as well.
#
# written by Aaron Lun
# a long time ago. 
# last modified 22 March 2017
{
    # Verifying patterns and overhangs.
    overhang <- rep(as.integer(overhang), length.out=length(pattern))
    all.patterns <- vector("list", length(pattern))
    remainder <- integer(length(overhang))

    for (p in seq_along(pattern)) {
        ps <- DNAString(pattern[p])
        if (reverseComplement(ps)!=ps) { 
            stop("recognition site must be an inverse palindrome") 
        }
        all.patterns[[p]] <- ps

        oh <- overhang[p]
        if (oh > length(ps) || oh < 0L) {
            stop("overhang must be a non-negative integer not greater than pattern length") 
        }

        even <- (oh %% 2L)==0L
        if (even != (length(ps)%%2L==0L)) {
            stop("both 'overhang' and pattern length must be either even or odd")
        }

        remainder[p] <- (length(ps) - oh)/2L
    }

    # Checking BSgenome object. 
    if (is(bs, "BSgenome")) {
        ref.names <- seqnames(bs)
        gen <- genome(bs)
    } else {
        bs <- readDNAStringSet(bs)
        ref.names <- names(bs)
        gen <- NA
    }

    # Running through the options.
    nchrs <- length(ref.names)
    original <- vector("list", nchrs)
    for (i in seq_len(nchrs)) {
        chr <- ref.names[i]
        cur.seq <- bs[[chr]]
        chrlen <- length(cur.seq)

        combined <- NULL
        for (p in seq_along(all.patterns)) {
            x <- matchPattern(all.patterns[[p]], cur.seq, fixed="subject")
            match.start <- start(x)
            starts <- match.start + remainder[p]
            ends <- match.start + remainder[p] - 1L + overhang[p]
            frags <- GRanges(chr, IRanges(c(1L, starts), c(ends, chrlen)))

            if (is.null(combined)) {
                combined <- frags
            } else {
                olap <- findOverlaps(combined, frags)
                combined <- pintersect(combined[queryHits(olap)], frags[subjectHits(olap)])
            }
        }

        combined$hit <- NULL
        original[[i]] <- combined
    }

    suppressWarnings(original <- do.call(c, original))
    seqlevels(original) <- ref.names
    suppressWarnings(seqinfo(original) <- Seqinfo(ref.names, seqlengths=seqlengths(bs), genome=gen))
    sort(original)
}
