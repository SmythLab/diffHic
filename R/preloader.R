preloader <- function(fnames, param=NULL, retain=NULL)
# This function produces a list of functions that allow data for each chromosome pair in each file to be loaded.
# The idea is to produce the functions and then to call them within a loop for each function as required.
{
    if (is.null(param)) { 
        chrs <- frag.by.chr <- discard <- restrict <- NULL
        do.restrict <- paired <- FALSE
        cap <- NA_integer_
    } else {
        # Checking how we want to do restriction.
        restrict <- param$restrict
    	do.restrict <- length(restrict)!=0L
    	paired <- attributes(restrict)$paired
    	if (is.null(paired)) { 
            paired <- FALSE 
        } else if (paired) { 
            nr <- length(restrict)/2L
        }
        
        # Pulling out the chromosomes we want.
        fragments <- param$fragments
        if (!.isDNaseC(fragments=fragments)) { 
            cap <- param$cap
            frag.by.chr <- .splitByChr(fragments) 
            chrs <- frag.by.chr$chrs
        } else {
            cap <- NA_integer_ # pointless to cap as all IDs are zeroes.
            all.lengths <- seqlengths(fragments)
            chrs <- names(all.lengths)
            first <- last <- integer(length(chrs)) # set to zero if we're not binning.
            names(first) <- names(last) <- chrs
            frag.by.chr <- list(first=first, last=last)
        }
    
        # Pulling out the discard ranges.
        discard <- .splitDiscards(param$discard)
    }

    # Iterating across all chromosome pairs and creating functions for extraction.
    y <- path.expand(fnames)
	overall <- list()
	for (ix in seq_along(y)) { 
		all.data <- loadChromos(y[ix])

        # Initial restriction on both entries being present in the list.
        if (do.restrict) {
            keep <- all.data$anchor1 %in% restrict & all.data$anchor2 %in% restrict
            all.data <- all.data[keep,]
        }
		current <- split(all.data$anchor2, all.data$anchor1)

        # Throwing an error if we see a chromosome not present in 'chrs'.
        if (!is.null(chrs)) { 
            leftovers <- setdiff(c(all.data$anchor1, all.data$anchor2), chrs)
            if (length(leftovers)) {
                stop("'", leftovers[1], "' in '", y[ix], "' is not present in fragment chromosomes")
            }
        }

		for (ac in names(current)) {
			subcurrent <- current[[ac]]

			# Paired restriction if necessary.
			if (do.restrict && paired) {
                tmatch <- which(restrict==ac)
                tmatch <- tmatch + nr * (-1L)^(tmatch > nr) # Getting the chromosome in the other column.
                subcurrent <- intersect(subcurrent, restrict[tmatch])
                if (!length(subcurrent)) { 
                    next 
                }
            }
			
            if (is.null(overall[[ac]])) { 
                overall[[ac]] <- list() 
            }
			for (tc in subcurrent) {
				if (is.null(overall[[ac]][[tc]])) { 
                    overall[[ac]][[tc]] <- rep(list(EMPTY_FUN(retain)), length(y)) 
                }
				overall[[ac]][[tc]][[ix]] <- RETRIEVAL_FUN(y[ix], ac, tc, frag.by.chr, discard, cap, retain=retain)
			}
		}
	}
	return(overall)
}

EMPTY_FUN <- function(retain) {
    if (!is.null(retain)) { 
        out <- matrix(0L, 0, length(retain)) 
        colnames(out) <- retain
        out <- data.frame(out)
        return(function() { out })
    } else {
        return(NULL)
    }
}

RETRIEVAL_FUN <- function(fname, anchor1, anchor2, chr.limits, discard, cap, retain) { 
    adisc <- discard[[anchor1]]
    tdisc <- discard[[anchor2]]

    check.limits <- !is.null(chr.limits)
    first1 <- chr.limits$first[[anchor1]]
    first2 <- chr.limits$first[[anchor2]]
    last1 <- chr.limits$last[[anchor1]] 
    last2 <- chr.limits$last[[anchor2]] 

    # Evaluating promises.
    force(fname)
    force(cap) 
    force(retain)

    function() {
        # Pulling out the reads, binning if necessary, and checking fidelity of the input.
        out <- .getPairs(fname, anchor1, anchor2)
        dim(out$anchor1.id) <- dim(out$anchor2.id) <- NULL
        check <- .Call(cxx_check_input, out$anchor1.id, out$anchor2.id)
    
        # Checking that we're all on the right chromosome.
        if (check.limits && nrow(out)) { 
            if (max(out$anchor1.id) > last1 || min(out$anchor1.id) < first1) {
                stop("anchor1 index outside range of fragment object") 
            }
            if (max(out$anchor2.id) > last2 || min(out$anchor2.id) < first2) {
                stop("anchor2 index outside range of fragment object") 
            }
        }
    
        # Overlapping with those in the discard intervals.
        if (!is.null(adisc) || !is.null(tdisc)) {
            a.hits <- t.hits <- FALSE
            if (!is.null(adisc)) {
                a.hits <- overlapsAny(IRanges(out$anchor1.pos, out$anchor1.pos+abs(out$anchor1.len)-1L), adisc, type="within")
            }
            if (!is.null(tdisc)) { 
                t.hits <- overlapsAny(IRanges(out$anchor2.pos, out$anchor2.pos+abs(out$anchor2.len)-1L), tdisc, type="within")
            }
            out <- out[!a.hits & !t.hits,,drop=FALSE]
        }
    
        # Removing read pairs above the cap for each restriction fragment pair.
        if (!is.na(cap)) { 
            capped <- .Call(cxx_cap_input, out$anchor1.id, out$anchor2.id, cap)
            out <- out[capped,]
        }
    
        if (!is.null(retain)) { 
            out <- out[,retain]
        } 
        return(out)
    }
}

