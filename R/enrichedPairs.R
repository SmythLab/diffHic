enrichedPairs <- function(data, flank=5, exclude=0, assay.in=1, assay.out=NULL)
# For each bin pair in 'data', this function counts the number of read pairs in 
# the neighbouring bin pairs in 'data', with four defined neighbourhood types.
# 
# written by Aaron Lun
# created 23 April 2014
# last modified 2 March 2017
{
	flank <- as.integer(flank)
	exclude <- as.integer(exclude)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	if (exclude < 0L) { stop("exclude width must be a positive integer") }
	if (flank <= exclude) { stop("exclude width must be less than the flank width") }
    .check_StrictGI(data)
	
	rdata <- .splitByChr(regions(data))
	last.id <- rdata$last
	first.id <- rdata$first
	
	# Running through each pair of chromosomes.
	np <- nrow(data)
    nl <- ncol(data)
	all.chrs <- as.character(seqnames(regions(data)))
	aid <- anchors(data, type="first", id=TRUE)
	by.chr <- split(seq_len(np), all.chrs[aid])
	tid <- anchors(data, type="second", id=TRUE)

    # Setting up the output for each mode.
    modes <- .neighbor_locales(assay.out)
    count.output <- lapply(modes, FUN=function(x) matrix(0L, nrow=np, ncol=nl))
    n.output <- lapply(modes, FUN=function(x) numeric(np))

	for (anchor in names(by.chr)) {
		next.chr <- by.chr[[anchor]]
		next.chr <- split(next.chr, all.chrs[tid[next.chr]])
		a.len <- last.id[[anchor]] - first.id[[anchor]] + 1L

		for (target in names(next.chr)) {
			current.pair <- next.chr[[target]]
			all.a <- aid[current.pair] - first.id[[anchor]] 
			all.t <- tid[current.pair] - first.id[[target]]
			t.len <- last.id[[target]] - first.id[[target]] + 1L
			
            o <- order(all.a, all.t)
            all.a <- all.a[o]
            all.t <- all.t[o]
            all.c <- assay(data, assay.in)[current.pair,,drop=FALSE][o,,drop=FALSE]

            # Getting counts for each library and type of neighbouring region.
            for (lib in seq_len(nl)) { 
                collected <- .Call(cxx_quadrant_bg, all.a, all.t, all.c[,lib],
                                   flank, exclude, a.len, t.len, anchor==target)

                for (m in seq_along(modes)) {
                    cur.counts <- collected[[1]][[m]]
                    cur.counts[o] <- cur.counts
                    count.output[[m]][current.pair,lib] <- cur.counts
                }

                if (lib==1L) {
                    for (m in seq_along(modes)) {
                        cur.n <- collected[[2]][[m]]
                        cur.n[o] <- cur.n
                        n.output[[m]][current.pair] <- cur.n
                    }
                }
            }
		}
	}

    # Adding the counts to the data and returning the object.
    n.names <- .neighbor_numbers(modes)
    for (m in seq_along(modes)) { 
        assay(data, modes[m]) <- count.output[[m]]
        mcols(data)[[n.names[m]]] <- n.output[[m]]
    }
	return(data)
}

.neighbor_locales <- function(x=NULL) { 
    if (is.null(x)) { 
        modes <- c("quadrant", "vertical", "horizontal", "surrounding") 
    } else { 
        modes <- x
        if (!is.character(modes) || length(modes)!=4L || anyDuplicated(modes)) { 
            stop("'assay.out' must be a character vector with 4 unique names")
        }
    }
    return(modes)
}

.neighbor_numbers <- function(x=NULL) { 
    if (is.null(x)) x <- .neighbor_locales()
    paste0("N.", x) 
}
