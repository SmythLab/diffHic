###################################################################################################
# This provides some functions to simulation fragments and the count directories.

simgen <- function(dir, num, chromos) {
	bonus<-c(0L, cumsum(as.integer(unlist(chromos))))
	overall<-NULL
	for (i in 1:length(chromos)) { 
		max.anchor<-chromos[[i]]
		for (j in 1:i) {
			max.target<-chromos[[j]]
			anchors<-as.integer(floor(runif(num, 1, max.anchor)))
			targets<-as.integer(floor(runif(num, 1, max.target)))
			if (i==j){
				anchor.1<-pmax(anchors, targets)
				target.1<-pmin(anchors, targets)
				anchors<-anchor.1
				targets<-target.1
			}
			cyrrebt<-data.frame(anchor.id=anchors+bonus[i], target.id=targets+bonus[j], 
					anchor.pos=0L, target.pos=0L, anchor.len=0L, target.len=0L)
			overall<-rbind(overall, cyrrebt)
		}
	}
	tmpfrags<-GRanges(rep(names(chromos), chromos), IRanges(1:sum(chromos), 1:sum(chromos)))
	savePairs(overall, dir, pairParam(tmpfrags))
}

# Spawning a new cut site set-up.

simcuts<-function(chromos, min=500, max=2000, overlap=0L) {
	cuts<-list()
	overlap <- as.integer(overlap)
	for (i in 1:length(chromos)) { 
		frags<-as.integer(runif(chromos[[i]], min, max)) 
		frag.ends<-cumsum(frags) - 0:(chromos[[i]]-1L)*overlap
		frag.starts<-c(1L, frag.ends[-chromos[[i]]]+1L-overlap)
		cur_chr<-names(chromos)[i]
		cuts[[cur_chr]]<-GRanges(cur_chr, IRanges(frag.starts, frag.ends))
		seqlengths(cuts[[cur_chr]]) <- max(frag.ends)
	}
	names(cuts)<-NULL
	suppressWarnings(cuts<-do.call(c, cuts))
	return(cuts)
}

# Adding in some positional information to each HDF5 file.

augmentsim <- function(infile, frags, rlen=10) {
	allfs <- start(frags)
	allfe <- end(frags)
   	x <- h5ls(infile)
	x <- x[x$otype=="H5I_DATASET",]

	everything <- list()
	for (i in 1:nrow(x)) { 
		cpath <- file.path(x$group[i], x$name[i])
		collected <- h5read(infile, cpath)
		num <- nrow(collected)
		
		a.s <- allfs[collected$anchor.id]		
		a.e <- allfe[collected$anchor.id]
		astr <- rbinom(num, 1, 0.5)==1L
		a.s2 <- ifelse(astr, a.s, pmin(a.e, a.s - rlen + 1L))
		a.e2 <- ifelse(astr, pmax(a.s, a.e - rlen + 1L), a.e)
		collected$anchor.pos <- as.integer(runif(nrow(collected), min=a.s2, max=a.e2))
		collected$anchor.len <- as.integer(rlen * ifelse(astr, 1, -1))

		t.s <- allfs[collected$target.id]
		t.e <- allfe[collected$target.id]
		tstr <- rbinom(num, 1, 0.5)==1L
		t.s2 <- ifelse(astr, t.s, pmin(t.e, t.s - rlen + 1L))
		t.e2 <- ifelse(astr, pmax(t.s, t.e - rlen + 1L), t.e)
		collected$target.pos <- as.integer(runif(nrow(collected), min=t.s2, max=t.e2))
		collected$target.len <- as.integer(rlen * ifelse(tstr, 1, -1))

		everything[[i]] <- collected
	}

	savePairs(do.call(rbind, everything), infile, pairParam(frags))
}

# Discard data.

makeDiscard <- function(ndisc, sizeof, chromosomes) {
	chosen <- sample(length(chromosomes), ndisc, replace=TRUE)
	chosen.pos <- runif(ndisc, 1, chromosomes[chosen]-sizeof)
	reduce(GRanges(names(chromosomes)[chosen], IRanges(chosen.pos, chosen.pos+sizeof-1L)))
}

###################################################################################################

