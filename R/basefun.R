.splitByChr <- function(ranges)
# Gets the start and end for each chromosome in the sorted GRanges. 
{
    chrs <- as.character(runValue(seqnames(ranges)))
    if (anyDuplicated(chrs)) { stop("ranges for each chromosome should be consecutive") }
    ref.len <- runLength(seqnames(ranges))
    end.index <- cumsum(ref.len)
    start.index <- end.index - ref.len + 1L
    names(end.index) <- names(start.index) <- chrs
    return(list(chrs=chrs, first=start.index, last=end.index))
}

.splitDiscards <- function(discard) 
# Splits the discard GRanges into a list of constituent chromosomes,
# along with IRanges for everything. This allows easy access to the
# ranges on individual chromosomes, rather than overlapping with everything.
{
    if (is.null(discard) || length(discard)==0L) { return(NULL) }
    discard <- sort(discard)
    all.chrs <- as.character(runValue(seqnames(discard)))
    all.len <- runLength(seqnames(discard))
    chr.ends <- cumsum(all.len)
    chr.starts <- c(1L, chr.ends[-length(chr.ends)]+1L)

    nchrs <- length(all.chrs)
    output <- vector("list", nchrs)
    for (i in seq_len(nchrs)) { 
        chr <- all.chrs[i]
        ix <- chr.starts[i]:chr.ends[i]
        output[[chr]] <- reduce(ranges(discard[ix]))
    }

    return(output)
}

####################################################################################################

.isDNaseC <- function(param, fragments) {
    if (missing(fragments)) { 
        fragments <- param$fragments
    }
    return(length(fragments)==0L)
}

.assignBins <- function(param, width, restricted=FALSE)
# Determines which bin each restriction fragment is in. Also records the rounded
# start and stop site for each bin. Returns a set of bin ids for each restriction
# fragment on each chromosome, as well as the coordinates of each bin.
{
    fragments <- param$fragments
    width <- as.integer(width)
    out.ids <- integer(length(fragments))
    last <- 0L 

    frag.data <- .splitByChr(fragments)
    out.ranges <- nfrags <- vector("list", length(frag.data$chrs))

    # Ignoring fragments not in our restricted set.
    if (restricted && length(param$restricted)) { 
        chrs.to.use <- which(frag.data$chrs %in% param$restricted)
    } else {
        chrs.to.use <- seq_along(frag.data$chrs)
    }

    for (x in chrs.to.use) { 
        curindex <- frag.data$first[x]:frag.data$last[x]
        curf <- fragments[curindex]
        mids <- (start(curf)+end(curf))/2
        bin.id <- as.integer((mids-0.1)/width)+1L 
        # The '-0.1' in the preceding step reduces 'mids' that are exact multiples 
        # of 'width', so each bin is from (n*width, (n+1)*width] for integer 'n'.

        # Setting up the bin ID for each fragment.
        processed <- rle(bin.id)
        ns <- length(processed$value)
        processed$values <- seq_len(ns)
        nfrags[[x]] <- processed$length
        out.ids[curindex] <- inverse.rle(processed)+last
        last <- last+ns
        
        # Setting up the restriction-fragment-rounded coordinates for each bin.
        endx <- cumsum(processed$length)
        startx <- rep(1L, ns)
        if (ns>=2L) { 
            startx[-1] <- endx[-ns]+1L 
        }
        out.ranges[[x]] <- GRanges(frag.data$chrs[x], IRanges(start(curf[startx]), end(curf[endx])))
    }

    # Wrapping up.
    suppressWarnings(out.ranges <- do.call(c, out.ranges))
    seqlevels(out.ranges) <- seqlevels(fragments)
    seqlengths(out.ranges) <- seqlengths(fragments)
    out.ranges$nfrags <- unlist(nfrags)
    return(list(id=out.ids, region=out.ranges))
}

.createBins <- function(param, width, restricted=FALSE) 
# This creates regular contiguous bins of size 'width'. Each bin
# is assigned to itself; allocation of read pairs into bins is done below.
# This allows free-floating bins for use with DNase-C data.
{
    ref.len <- seqlengths(param$fragments)
    chr.names <- names(ref.len)
    everything <- vector("list", length(ref.len))

    # Ignoring fragments not in our restricted set.
    if (restricted && length(param$restricted)) { 
        chrs.to.use <- which(chr.names %in% param$restricted)
    } else {
        chrs.to.use <- seq_along(chr.names)
    }

    for (i in chrs.to.use) {
        chr.len <- ref.len[i]
        bin.dex <- seq_len(ceiling(chr.len/width))
        end.pt <- pmin(bin.dex * width, chr.len)
        current <- GRanges(chr.names[i], IRanges((bin.dex - 1L)*width + 1L, end.pt))
        everything[[i]] <- current
    }

    suppressWarnings(everything <- do.call(c, everything))
    seqlengths(everything) <- ref.len
    everything$nfrags <- 0L
    return(list(id=seq_along(everything), region=everything))
}

####################################################################################################

.binReads <- function(pairs, width, first1, first2, last1, last2)
# Binning the read pairs into bins of size 'width',
# based on the 5' coordinates of each read.
{
    anchor1.id <- .readToBin(pairs$anchor1.pos, pairs$anchor1.len, 
                             bin.width=width, first.bin=first1, last.bin=last1)
    anchor2.id <- .readToBin(pairs$anchor2.pos, pairs$anchor2.len, 
                             bin.width=width, first.bin=first2, last.bin=last2)

    new.pairs <- data.frame(anchor1.id=anchor1.id, anchor2.id=anchor2.id)
    new.pairs <- .enforcePairOrder(new.pairs)
    o <- order(new.pairs$anchor1.id, new.pairs$anchor2.id)
    new.pairs[o,]
}

.readToBin <- function(read.pos, read.len, bin.width, first.bin, last.bin)
# Converts read positions to bin IDs. If any reads have start positions in 
# bins beyond the last.bin, these are preserved to ensure they get caught
# by the error-checking machinery later on.
{
    id <- as.integer(ceiling(read.pos/bin.width)) - 1L + first.bin
    is.ok.rev <- read.len < 0L & id <= last.bin # i.e., reverse reads with valid starting positions.

    rev5 <- read.pos[is.ok.rev] - read.len[is.ok.rev] - 1L
    rid <- as.integer(ceiling(rev5/bin.width)) - 1L + first.bin
    id[is.ok.rev] <- pmin(last.bin, rid) # pmin only rescues reverse reads that were already okay.   
    return(id)
}

.enforcePairOrder <- function(pairs) 
# Enforcing 1 >= 2.
{
    swap <- pairs$anchor2.id > pairs$anchor1.id
    if (any(swap)) { 
        flipping <- pairs[swap,]
        has.one <- grepl("1", colnames(flipping))
        has.two <- grepl("2", colnames(flipping))
        colnames(flipping)[has.one] <- sub("1", "2", colnames(flipping)[has.one])
        colnames(flipping)[has.two] <- sub("2", "1", colnames(flipping)[has.two])
        for (field in colnames(flipping)) { 
            pairs[[field]][swap] <- flipping[[field]]
        }
    }
    return(pairs)
}

####################################################################################################

