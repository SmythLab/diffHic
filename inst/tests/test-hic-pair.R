###################################################################################################
# This script is designed to test the pair-identifying capabilities of the hiC machinery i.e. preparePairs. 
# We start with unit tests for individual components of the preparePairs C++ code.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

# Checking CIGAR.

checkCIGAR <- function(cigar, rstrand) {
	out <- .Call(diffHic:::cxx_test_parse_cigar, cigar, rstrand)
	if (is.character(out)) { stop(out) }
	
	true.alen <- GenomicAlignments::cigarWidthAlongReferenceSpace(cigar)
	if (out[1]!=true.alen) { stop("mismatch in alignment length") }

	offset <- 0
	as.rle <- GenomicAlignments::cigarToRleList(cigar)[[1]]
	if (rstrand) { as.rle <- rev(as.rle) }
	if (runValue(as.rle)[1]=="H") { offset <- runLength(as.rle)[1] } 
	if (out[2]!=offset) { stop("mismatch in offsets") }

	return(c(alen=true.alen, offset=offset))
}

checkCIGAR("5H20M", TRUE)
checkCIGAR("5H20M", FALSE)

checkCIGAR("5H20M6H", TRUE)
checkCIGAR("5H20M6H", FALSE)

checkCIGAR("5H20M6S", TRUE)
checkCIGAR("5H20M6S", FALSE)

checkCIGAR("20M5I30M", TRUE)
checkCIGAR("20M5D30M", TRUE)
checkCIGAR("20M5N30M", TRUE)
checkCIGAR("20M5P30M", TRUE)

checkCIGAR("20M5X30M", TRUE)
checkCIGAR("20M5=30M", TRUE)

checkCIGAR("10M2I5M3D20M", TRUE)
checkCIGAR("10M2I3D5N20M", TRUE)
checkCIGAR("10M2P5M3D20I", TRUE)
checkCIGAR("10M2I3D5N20D", TRUE)

checkCIGAR("5H11M5D8I9X6=6H", TRUE)
checkCIGAR("5H20X5D8P9N72M", FALSE)
checkCIGAR("5S19M34=55D8X20M6S", TRUE)
checkCIGAR("5S1M3=5D18X2M6S", FALSE)

# Checking fragment assignment.

assign2fragment <- function(starts, ends, chr, pos, rstrand, len) {
	out <- .Call(diffHic:::cxx_test_fragment_assign, starts, ends, chr, pos, rstrand, len)
	if (is.character(out)) { stop(out) }

	chr <- chr + 1L
	if (rstrand) { 
		fiveprime <- min(pos + len -1L, tail(ends[[chr]], 1))
		stopifnot(ends[[chr]][out] >= fiveprime && (out==1L || ends[[chr]][out-1] < fiveprime))
	} else { 
		fiveprime <- pos 
		stopifnot(starts[[chr]][out] <= fiveprime && (out==length(starts[[chr]]) || starts[[chr]][out+1] > fiveprime))
	}
	
	out
}

starts <- list( c(1L, 100L, 200L, 300L, 400L, 500L), # chr1
                c(1L, 100L, 200L, 300L, 400L, 500L)) # chr2
ends <- list( c(103L, 203L, 303L, 403L, 503L, 1000L), # chr1
              c(103L, 203L, 303L, 403L, 503L, 1000L)) # chr2

assign2fragment(starts, ends, 0L, 94L, TRUE, 10L)
assign2fragment(starts, ends, 0L, 95L, TRUE, 10L)
assign2fragment(starts, ends, 0L, 99L, FALSE, 10L)
assign2fragment(starts, ends, 0L, 100L, FALSE, 10L)

assign2fragment(starts, ends, 0L, 203L, TRUE, 1L)
assign2fragment(starts, ends, 0L, 204L, TRUE, 1L)
assign2fragment(starts, ends, 0L, 209L, FALSE, 1L)
assign2fragment(starts, ends, 0L, 300L, FALSE, 1L)

assign2fragment(starts, ends, 1L, 1L, FALSE, 10L)
assign2fragment(starts, ends, 1L, 991L, FALSE, 10L)
assign2fragment(starts, ends, 1L, 992L, FALSE, 10L)
assign2fragment(starts, ends, 1L, 991L, TRUE, 10L)
assign2fragment(starts, ends, 1L, 992L, TRUE, 10L)

starts <- list( c(1L, 1L, 100L, 200L, 300L, 400L, 500L, 997L), # chr1, with nesting at the start and end.
                c(1L, 1L, 100L, 200L, 300L, 400L, 500L, 997L)) # chr2
ends <- list( c(4L, 103L, 203L, 303L, 403L, 503L, 1000L, 1000L), # chr1
              c(4L, 103L, 203L, 303L, 403L, 503L, 1000L, 1000L)) # chr2

assign2fragment(starts, ends, 0L, 1L, FALSE, 10L)
assign2fragment(starts, ends, 0L, 1L, TRUE, 10L)
assign2fragment(starts, ends, 0L, 991L, FALSE, 10L)
assign2fragment(starts, ends, 0L, 991L, TRUE, 10L)
try(assign2fragment(starts, ends, 0L, 1000L, TRUE, 10L)) # This should fail, as it gets assigned into the nested fragment.

###################################################################################################
# We also set up a full simulation for the entire function.

suppressPackageStartupMessages(require("rhdf5"))
source("simsam.R")

comp<-function (fname, npairs, max.cuts, sizes=c(100, 500), singles=0, rlen=10, spacer=rlen, 
		yield=max(1L, round(npairs/runif(1, 2, 10)))) {
	rlen<-as.integer(rlen)
	spacer<-as.integer(spacer)
	if (min(sizes) <= rlen) { stop("min fragment must be greater than read length") } 
	# Necessary for proper assignment, especially at the start of the chromosome when reverse 
	# reads are bounded at zero (i.e. their 5' ends would not be defined if 1+rlen > fragmentsize)

	# Randomly generating fragment lengths for the chromosome.
	fragments<-list()
	cut.starts<-list()
	chromosomes<-rep(0L, length(max.cuts));
	outfrags<-list()
	for (i in 1:length(max.cuts)) { 
		fragments[[i]]<-as.integer(round(runif(max.cuts[[i]], sizes[1], sizes[2])))
		ends<-cumsum(fragments[[i]])+1:max.cuts[[i]]*spacer
		if (max.cuts[[i]]>1) { 
			cut.starts[[i]]<-c(1L, (ends+1L)[1:(max.cuts[[i]]-1)])
		} else {
			cut.starts[[i]]<-1L
		}
		outfrags[[i]]<-GRanges(names(max.cuts)[i], IRanges(cut.starts[[i]], ends-spacer))
		chromosomes[i]<-tail(ends, 1)-spacer;
	}
	suppressWarnings(outfrags<-do.call(c, outfrags))
	names(chromosomes)<-names(max.cuts)
	names(fragments)<-names(max.cuts)
	names(cut.starts)<-names(max.cuts)

	# Randomly generating reads (a la getPETSizes' example).
    names<-paste('x', rep(1:npairs, 2), sep=".");
    chrs<-sample(length(chromosomes), length(names), replace=TRUE);
    pos<-integer(length(names));
	frag.ids<-integer(length(names))
    str<-rbinom(length(names), 1, 0.5)==1; 

    # Assigning positions to all of them. Some finesse is necessary when 
	# there is no spacer (i.e. original genome with positive overhang). We
	# still allow reads to span restriction sites, though.
    for (i in 1:length(chromosomes)) {
        current<-chrs==i;
		chosen.frags<-as.integer(runif(sum(current), 1, length(fragments[[i]]))+0.5)
		frag.ids[current]<-chosen.frags
		my.ends <- cut.starts[[i]] + fragments[[i]]
		if (spacer > 0) {
			pos[current]<-as.integer(runif(sum(current), cut.starts[[i]][chosen.frags], my.ends[chosen.frags]-rlen))
		} else {
			cur.for <- str[current]
			forward.frag <- chosen.frags[cur.for]
			pos[current][cur.for] <- as.integer(runif(sum(cur.for), 
					cut.starts[[i]][forward.frag],
					ifelse(forward.frag!=length(cut.starts[[i]]),
						cut.starts[[i]][forward.frag+1L],
						my.ends[forward.frag])))
			
			reverse.frag <- chosen.frags[!cur.for]
			starter <- integer(sum(!cur.for))
			possible.zero <- reverse.frag==1L
			starter[possible.zero] <- 1L
			starter[!possible.zero] <- my.ends[reverse.frag[!possible.zero]-1L] 
			pos[current][!cur.for] <- pmax(1L, 
				as.integer(runif(sum(!cur.for), starter, my.ends[chosen.frags[!cur.for]])) - rlen + 1L)
			# The boundary conditions here are tricky, but because runif will never generate values
			# equal to the boundary, rounding down will always ensure a position in the correct fragment.
		}
    }

    # Throwing them into the SAM file generator. 
	reversi<-c(1:npairs+npairs, 1:npairs)
    out<-simsam(fname, names(chromosomes)[chrs], pos, str, chromosomes, names=names, len=rlen, 
			is.first=c(rep(TRUE, npairs), rep(FALSE, npairs)), is.paired=TRUE,
			mate.chr=names(chromosomes)[chrs][reversi], mate.pos=pos[reversi], mate.str=str[reversi])

	if (singles) { 
	    # Adding some singles. You'll get some warnings regarding spacer regions as we're not finessing it.
    	snames<-schrs<-spos<-NULL
        snames<-paste('y', 1:singles, sep=".");
        schrs<-sample(length(chromosomes), singles, replace=TRUE);
        spos<-integer(singles);
		for (i in 1:length(chromosomes)) {
	       	scurrent<-schrs==i;
        	spos[scurrent]<-as.integer(round(runif(sum(scurrent), 1, chromosomes[i])))
		}
	
    	tempname<-file.path(dir, "temp")
    	sstr<-rbinom(singles, 1, 0.5)==1; 
		out2<-simsam(tempname, names(chromosomes)[schrs], spos, sstr, chromosomes, names=snames, len=rlen)
		more.temp<-file.path(dir, "temp2")
		out<-mergeBam(c(out, out2), more.temp, indexDestination=TRUE, overwrite=TRUE)
		file.rename(more.temp, out)
	}

	# Resorting by name.
	temp<-sortBam(out, "temp", byQname=TRUE)
	file.rename(temp, out)

	################ THEORETICAL MATCH ###################

	# This gets the status: 0 for okay, 1 for dangling end, 3 for self.circles and and 2 for other stuff.
	getstatus <- function(chr1, pos1, frag1, str1, chr2, pos2, frag2, str2) {
		codes<-rep(0L, length(chr1));
		potentials<-chr1==chr2 & frag1==frag2
		
		same.str<-str1==str2;
		codes[ same.str & potentials ] <- 2L;
		potentials <- potentials & !same.str;
		
		self.circle<-(str1 & pos2+rlen <= pos1) | (str2 & pos1+rlen <= pos2);
		codes[ self.circle & potentials ] <- 3L
		potentials<-potentials & !self.circle
		
		overextension<-(str1 & pos1 > pos2) | (str2 & pos2 > pos1) # Don't need to check +rlen, as they're al the same.
		codes[ overextension & potentials ]<- 2L;
		potentials<-potentials & !overextension;
		
		codes[potentials] <- 1L;
		return(codes);
	}

	# This gets the distance from the 5' end of the read to the next restriction site that it is pointing to (past the one it covers, if it is incomplete).
	getlen <- function(chr, pos, start, str) {
		len.out<-rep(0L, length(chr));
		for (x in 1:length(cut.starts)) { 
			chosen<-x==chr;
			dist2cut<-ifelse(str[chosen], cut.starts[[x]][start[chosen]]+fragments[[x]][start[chosen]]-pos[chosen],
					pos[chosen]+rlen-cut.starts[[x]][start[chosen]])
			len.out[chosen]<-dist2cut;
		}
		return(len.out);
	}	

	# Now, actually assembling the theoretical values.
	primary<-1:npairs
	secondary<-npairs+1:npairs
	pchrs<-chrs[primary]
	ppos<-pos[primary]
	pfrag<-frag.ids[primary]
	pstr<-str[primary]
	schrs<-chrs[secondary]
	spos<-pos[secondary]
	sfrag<-frag.ids[secondary]
	sstr<-str[secondary]
	codes<-getstatus(pchrs, ppos, pfrag, pstr, schrs, spos, sfrag, sstr)
	frag.lens<-getlen(pchrs, ppos, pfrag, pstr)+getlen(schrs, spos, sfrag, sstr)
	panchor <- pchrs > schrs | 
		(pchrs==schrs & pfrag > sfrag) | 
		(pchrs==schrs & pfrag==sfrag & ppos > spos)
	inserts<-ifelse(pchrs==schrs, pmax(ppos, spos)-pmin(spos, ppos)+rlen, NA)
	orientations<-ifelse(pstr, 0L, ifelse(panchor, 1L, 2L))+ifelse(sstr, 0L, ifelse(panchor, 2L, 1L))

	################ ACTUAL MATCH ###################
	# Assembling the output list for comparison.

	tmpdir<-paste0(fname, "_temp")
	param <- pairParam(fragments=outfrags)
	diagnostics <- preparePairs(out, param, tmpdir, yield=yield);
   	
	stopifnot(sum(codes==1L)==diagnostics$same.id[["dangling"]])
	stopifnot(sum(codes==3L)==diagnostics$same.id[["self.circle"]])
	stopifnot(length(codes)==diagnostics$pairs[["total"]])
	stopifnot(singles==diagnostics$singles)

	# No support for testing chimeras, we use a fixed example below.
	stopifnot(diagnostics$unmapped.chimeras==0L) 
	stopifnot(diagnostics$chimeras[["total"]]==0L)
	stopifnot(diagnostics$chimeras[["mapped"]]==0L)
	stopifnot(diagnostics$chimeras[["invalid"]]==0L)

	# Anchor/target synchronisation is determined by order in 'fragments' (and thusly, in max.cuts).
	offset<-c(0L, cumsum(max.cuts))
	names(offset)<-NULL
	indices<-diffHic:::.loadIndices(tmpdir, seqlevels(outfrags))
	used<-indices
	fchrs<-as.character(seqnames(outfrags))

	for (i in 1:length(max.cuts)) {
		for (j in 1:i) {
			stuff<-(chrs[primary]==i & chrs[secondary]==j) | (chrs[primary]==j & chrs[secondary]==i) 
			stuff<-stuff & (codes==0L | codes==2L)
			pids<-frag.ids[primary][stuff];
			sids<-frag.ids[secondary][stuff];
			if (i > j) {
				which.is.which<-chrs[primary][stuff]==i & chrs[secondary][stuff]==j
 			    anchor<-ifelse(which.is.which, pids, sids)
				target<-ifelse(which.is.which, sids, pids)
			} else {
				anchor<-pmax(pids, sids)
				target<-pmin(pids, sids)
			}
			anchor<-anchor+offset[i]
			target<-target+offset[j]
			totes<-frag.lens[stuff]
			cur.ori<-orientations[stuff]
			cur.insert<-inserts[stuff]
			o<-order(anchor, target, totes, cur.ori, cur.insert)

			# Checking anchor/target/length/orientation/insert statistics (sorting to ensure comparability).
			achr<-names(max.cuts)[i]
			tchr<-names(max.cuts)[j]
			if (!(achr%in%names(indices)) || !(tchr %in% names(indices[[achr]]))) { 
				if (length(o)) { stop("true interactions are missing") }
				next
			}
			current<-h5read(tmpdir, file.path(achr, tchr))
			for (x in 1:ncol(current)) { attributes(current[,x]) <- NULL }
			collated <- diffHic:::.getStats(current, achr==tchr, outfrags)
			
			o2 <- order(current$anchor.id, current$target.id, collated$length, collated$orientation, collated$insert)
			stopifnot(identical(current$anchor.id[o2], anchor[o]))
			stopifnot(identical(current$target.id[o2], target[o]))
			stopifnot(identical(collated$length[o2], totes[o]))
			stopifnot(identical(collated$orientation[o2], cur.ori[o]))
			if (!identical(collated$insert[o2], cur.insert[o])) { 
				print(cbind(collated$insert[o2], cur.insert[o]))
				print(current[o2,])
			}
			stopifnot(identical(collated$insert[o2], cur.insert[o]))
				
			# Checking that we're looking at the right combination.
	        uniq.a<-unique(fchrs[current[,1]])
			uniq.t<-unique(fchrs[current[,2]])
			if (length(uniq.a)!=1L || length(uniq.t)!=1L) { stop("file contains more than one combination") }
			if (achr!=uniq.a || tchr!=uniq.t) { stop("file contains the incorrect combination") }
			used[[achr]][[tchr]]<-NULL
		}
	}

	# Checking there's nothing left.
	if (!is.null(unlist(used))) { stop("objects left unused in the directory") }

	# Length insert and orientation checking.
	keepers<-codes==0L | codes==2L
	valid.len<-frag.lens[keepers]
	valid.insert<-inserts[keepers]
	valid.ori<-orientations[keepers]
	o<-order(valid.len, valid.ori, valid.insert)
	auxiliary<-getPairData(tmpdir, param)
	o2<-do.call(order, auxiliary)
	if (!identical(valid.len[o], auxiliary$length[o2])) { stop("extracted fragment sizes don't match up") }
	if (!identical(valid.insert[o], auxiliary$insert[o2])) { stop("extracted inserts don't match up") }
	if (!identical(valid.ori[o], auxiliary$orientation[o2])) { stop("extracted orientations don't match up") }

	curdex <- h5ls(tmpdir)
	curdex <- curdex[curdex$otype=="H5I_DATASET",][1,]
	returned <- h5read(tmpdir, file.path(curdex$group, curdex$name))
	processed <- diffHic:::.getStats(returned, basename(curdex$group)==curdex$name, outfrags)
	return(head(data.frame(anchor.id=returned$anchor.id, target.id=returned$target.id, length=processed$length,
		orientation=processed$orientation, insert=processed$insert)))
}

####################################################################################################
# Initiating testing with something fairly benign.

set.seed(0)
dir<-"hic-test"
dir.create(dir)
fname<-file.path(dir, "out");
max.cuts<-c(chrA=20L, chrB=10L, chrC=5L)

comp(fname, npairs=20, max.cuts=max.cuts)
comp(fname, npairs=50, max.cuts=max.cuts)
comp(fname, npairs=100, max.cuts=max.cuts)

# Ramping up the aggression in terms of self-circles.

comp(fname, npairs=20, max.cuts=c(chrA=3L))
comp(fname, npairs=50, max.cuts=c(chrA=3L))
comp(fname, npairs=100, max.cuts=c(chrA=3L))

# Increasing the number of reads all round.

comp(fname, npairs=200, size=c(500, 1000), max.cuts=max.cuts);
comp(fname, npairs=1000, size=c(500, 1000), max.cuts=max.cuts);

# Adding some singletons.

comp(fname, npairs=200, size=c(20, 100), max.cuts=max.cuts, singles=10)
comp(fname, npairs=200, size=c(20, 100), max.cuts=max.cuts, singles=50)

# Making the fragments smaller.

comp(fname, npairs=200, size=c(60, 100), rlen=50, max.cuts=max.cuts)
comp(fname, npairs=500, size=c(60, 100), rlen=50, max.cuts=max.cuts)

# Trying out negative spacings (i.e., non-filled genomes)

comp(fname, npairs=200, size=c(500, 1000), spacer=-2, max.cuts=max.cuts);
comp(fname, npairs=1000, size=c(500, 1000), spacer=-2, max.cuts=max.cuts);
comp(fname, npairs=200, size=c(500, 1000), spacer=-4, max.cuts=max.cuts);
comp(fname, npairs=1000, size=c(500, 1000), spacer=-4, max.cuts=max.cuts);

comp(fname, npairs=200, size=c(20, 100), spacer=-2, max.cuts=max.cuts, singles=10)
comp(fname, npairs=200, size=c(20, 100), spacer=-2, max.cuts=max.cuts, singles=50)
comp(fname, npairs=200, size=c(20, 100), spacer=-4, max.cuts=max.cuts, singles=10)
comp(fname, npairs=200, size=c(20, 100), spacer=-4, max.cuts=max.cuts, singles=50)

comp(fname, npairs=200, size=c(60, 100), spacer=-4, rlen=50, max.cuts=max.cuts)
comp(fname, npairs=500, size=c(60, 100), spacer=-4, rlen=50, max.cuts=max.cuts)
comp(fname, npairs=200, size=c(60, 100), spacer=-2, rlen=50, max.cuts=max.cuts)
comp(fname, npairs=500, size=c(60, 100), spacer=-2, rlen=50, max.cuts=max.cuts)

# Trying it out with some more elements in a more restricted space.

comp(fname, npairs=500, size=c(100, 500), max.cuts=c(chrA=2L))
comp(fname, npairs=200, size=c(500, 1000), max.cuts=c(chrA=2L, chrB=1L))
comp(fname, npairs=1000, size=c(20, 50), max.cuts=c(chrA=1L))

# Adding lots of chromosomes.

max.cuts<-c(chrA=5L, chrB=6L, chrC=7L, chrD=4L, chrE=3L, chrF=1L, chrG=2L)
comp(fname, npairs=200, size=c(500, 1000), max.cuts=max.cuts);
comp(fname, npairs=1000, size=c(50, 100), max.cuts=max.cuts);
comp(fname, npairs=200, size=c(100, 500), max.cuts=max.cuts);
comp(fname, npairs=1000, size=c(200, 300), max.cuts=max.cuts);

###################################################################################################
# Trying to do simulations with chimeras is hellishly complicated, so we're just going to settle for 
# consideration of chimeras with a fixed example.

hic.file<-system.file("exdata", "hic_sort.bam", package="diffHic")
break.file<-system.file("exdata", "cuts.rds", package="diffHic")
cuts<-readRDS(break.file)
param <- pairParam(cuts)
tmpdir<-file.path(dir, "gunk")
cntdir<-file.path(dir, "gunkcount")

# Check out generator.R in inst/exdata for cross-referencing to names. 
# In order of appearance in printfun, all mapped read pairs are (for
# each data.frame in the list):

named <- list(c("chimeric.invalid.5", "good.1", "good.2", "good.3", "chimeric.invalid.4", "good.4", "chimeric.invalid.6", "other.1"),
		c("good.8", "good.5", "chimeric.good.4", "chimeric.good.5", "chimeric.good.1", "chimeric.invalid.2", "chimeric.invalid.1", "chimeric.invalid.3", "chimeric.invalid.7", "chimeric.good.2"),
		c("good.7", "good.6", "chimeric.good.3", "other.2"))

# We also have 3 unmapped reads, 3 dangling ends, 2 self-circles, 2 singletons.
# For chimeras, all have mapped 5' and 3' ends, 7 of which are invalid.
	
preparePairs(hic.file, param, tmpdir, dedup=FALSE)

printfun<-function(dir, named=NULL) {
	output<-list()
	ix <- 1L
	indices <- suppressWarnings(diffHic:::.loadIndices(tmpdir))
	for (ax in names(indices)) {
		if (is.null(output[[ax]])) { output[[ax]]<-list() }
		for (tx in names(indices[[ax]])) {
			extracted <- h5read(dir, file.path(ax, tx))
			processed <- diffHic:::.getStats(extracted, ax==tx, cuts)
			output[[ax]][[tx]] <- data.frame(anchor.id=extracted$anchor.id, target.id=extracted$target.id,
				length=processed$length, orientation=processed$orientation, insert=processed$insert)
			if (!is.null(named[[ix]])) { rownames(output[[ax]][[tx]])<-named[[ix]] }
			ix <- ix + 1L
		}
	}
	return(output)
}
printfun(tmpdir, named=named)

# Alright, so once duplicates are removed, we lose:
#  	self.1 (mapped/self.circles -> marked)
#   other.1 (mapped/other -> marked)	
#   dangling.1 (mapped/dangling -> marked)
#   chimeric.good.2 (mapped/chimeric$mapped/multi -> marked)
#	chimeric.invalid.1 (mapped/chimeric$mapped/multi/invalid -> marked)
#   chimeric.invalid.4 (--chimeric$multi/invalid)
#	chimeric.good.4 (mapped/chimeric$mapped/multi -> marked)
#	chimeric.good.5 (mapped/chimeric$mapped/multi -> marked)
#
# So, a gain of +7 to marked, a loss of -7 to mapped, a loss of -1 to each of
# dangling and self.circles, a loss of -4 to chimeric$mapped, -5 to
# chimeric$multi and -2 to chimeric$invalid. 
#
# Note that this will be considered the reference to which all downstream
# tests are compared, as dedup=FALSE is the default setting.

tmpdir2<-file.path(dir, "gunk2")
preparePairs(hic.file, param, tmpdir2)
named <- list(c("chimeric.invalid.5", "good.1", "good.2", "good.3", "chimeric.invalid.4", "good.4", "chimeric.invalid.6"),
	c("good.8", "good.5", "chimeric.good.1", "chimeric.invalid.2", "chimeric.invalid.3", "chimeric.invalid.7"),
	c("good.7", "good.6", "chimeric.good.3", "other.2"))
printfun(tmpdir2, named=named)

# Once invalid chimeras are removed, we see a loss of all rows corresponding to
# invalid chimeras in printfun. No change in the statistics should be observed.
# Note that chimeric.invalid.4 is still okay as ther invalid component is removed
# by duplicate removal (for some reason; that shouldn't happen in real data).

preparePairs(hic.file, param, tmpdir2, ichim=FALSE)
named <- list(c("good.1", "good.2", "good.3", "chimeric.invalid.4", "good.4"),
	c("good.8", "good.5", "chimeric.good.1"),
	c("good.7", "good.6", "chimeric.good.3", "other.2"))
printfun(tmpdir2, named=named)

# Throwing out those with poor mapping quality. We lose:
#   good.1 (mapped -> filtered)
#   good.3 (mapped -> filtered)
#   dangling.2 (mapped/dangling -> filtered)
#   dangling.3 (mapped/dangling -> filtered)
#   chimeric.good.1 (--chimeric$multi)
#	chimeric.invalid.2 (mapped/chimeric$mapped/multi/invalid -> filtered)
#   chimeric.invalid.3 (--chimeric$multi/invalid)
# 	chimeric.invalid.6 (mapped/chimeric$mapped/multi/invalid -> filtered)
#
# So, a gain of +6 for filtered, a loss of -6 for mapped, a loss of -2
# for dangling, a loss of -2 for chimeric$mapped, -4 for chimeric$multi
# and -3 for chimeric$invalid. For printfun, we see:
# 	A/A = chimeric.invalid.5, good.2, chimeric.invalid.4, good.4
#	B/A = good.8, good.5, chimeric.invalid.7
# 	B/B = good.7, good.6, chimeric.good.

preparePairs(hic.file, param, tmpdir2, minq=100)
named <- list(c("chimeric.invalid.5", "good.2", "chimeric.invalid.4", "good.4"),
	c("good.8", "good.5", "chimeric.good.1", "chimeric.invalid.3", "chimeric.invalid.7"),
	c("good.7", "good.6", "chimeric.good.3", "other.2"))
printfun(tmpdir2, named=named)

# Defining invalid chimeras based on distance instead of fragment ID. chimeric.good.1
# becomes an invalid chimera, as the distance between the 3' segment and the mate is
# 30 bp. This results in a +1 increase for invalid.chim.

preparePairs(hic.file, param, tmpdir2, chim.dist=20, ichim=FALSE)
named <- list(c("good.1", "good.2", "good.3", "chimeric.invalid.4", "good.4"),
	c("good.8", "good.5"),
	c("good.7", "good.6", "chimeric.good.3", "other.2"))
printfun(tmpdir2, named=named)

# chimeric.invalid.6 now becomes a valid chimera, as each pair of 5' end and 3' mate end
# is now a proper pair (inward-facing and less than chim.dist). -1 for invalid.chim.

preparePairs(hic.file, param, tmpdir2, chim.dist=2000, ichim=FALSE)
named <- list(c("good.1", "good.2", "good.3", "chimeric.invalid.4", "good.4", "chimeric.invalid.6"),
	c("good.8", "good.5", "chimeric.good.1"),
	c("good.7", "good.6", "chimeric.good.3", "other.2"))
printfun(tmpdir2, named=named)

###################################################################################################
# This tests whether they are counted properly. 

unlink(dir, recursive=TRUE) # Cleaning up

###################################################################################################

