normalizeCNV <- function(data, margins, prior.count=3, span=0.3, maxk=500, ...)
# This performs two-dimensional loess smoothing, using the counts and the 
# marginal counts to compute the abundance and the marginal fold-changes,
# respectively. Both are used as covariates in the model to smooth out any
# systematic differences in interaction intensity. The aim is to get rid
# of any CNV-induced bias, quantified by the differences in the marginals.
#
# written by Aaron Lun
# created 11 September 2014
# last modified 22 July 2015
{
	cont.cor <- 0.5
	cont.cor.scaled <- cont.cor * data$totals/mean(data$totals)
	ab <- aveLogCPM(counts(data), lib.size=data$totals, prior.count=cont.cor)
	mave <- aveLogCPM(counts(margins), lib.size=margins$totals, prior.count=prior.count)
	if (!identical(margins$totals, data$totals)) { 
		warning("library sizes should be identical for margin and data objects")
	}

	# Generating covariates.
	mab <- cpm(counts(margins), lib.size=margins$totals, log=TRUE, prior.count=prior.count) - mave
	matched <- matchMargins(data, margins)	
	ma.adjc <- mab[matched$amatch,,drop=FALSE] 
	mt.adjc <- mab[matched$tmatch,,drop=FALSE]

	offsets <- matrix(0, nrow=nrow(data), ncol=ncol(data))
	for (lib in seq_len(ncol(data))) {
		ma.fc <- ma.adjc[,lib]
		mt.fc <- mt.adjc[,lib]

		# Anchor/target distinction is arbitrary, so this coerces otherwise-identical 
		# points into the same part of the covariate space (see comment below).
		mfc1 <- (ma.fc + mt.fc)/2
		mfc2 <- abs(ma.fc - mt.fc)
		all.cov <- list(mfc1, mfc2, ab)
	
		# Fitting a loess surface with the specified covariates.	
		i.fc <- log2(counts(data)[,lib] + cont.cor.scaled[lib]) - ab 
		cov.fun <- do.call(lp, c(all.cov, nn=span, deg=1))
		fit <- locfit(i.fc ~ cov.fun, maxk=maxk, ..., lfproc=locfit.robust) 
		offsets[,lib] <- fitted(fit)
	}

	offsets <- offsets/log2(exp(1))
	offsets <- offsets - rowMeans(offsets)
	return(offsets)
}

##################### COMMENT ##########################
# You get two copy number changes for each bin pair, one for each region. The
# simplest coercion involves defining one covariate as the maximum change, and
# the other covariate as the minimum change. However, this gives a covariate
# space that is cut off past the diagonal. This won't be happily fitted in
# locfit, as it uses a rectangular grid (check out Computational Methods in
# Loader's book). You need all corners of the grid to interpolate, but one of
# those corners will be useless if it hangs on the wrong side of the diagonal.
# Instead, we rotate the space by 45 degrees, such that the diagonal is now a
# vertical line. The boundary of the grid now coincides with the true boundary
# of the space. All corners will now have sensible evaluations, such that 
# interpolation will be more reliable.
########################################################

matchMargins <- function(data, margins) 
# This function just matches the bin pairs in 'data' to the two indices of
# 'margins' that each bin corresponds to.
#
# written by Aaron Lun
# created 17 September 2014	
# last modified 20 March 2015 
{
	# Checking to ensure that the regions are the same.
	if (!identical(regions(data), regions(margins))) {
		stop("regions must be the same for bin pair and marginal counts") 
	}
	all.indices <- integer(length(regions(data)))
	id <- anchors(margins, id=TRUE)
	all.indices[id] <- seq_along(id)
	amatch <- all.indices[anchors(data, id=TRUE)]
	if (any(amatch==0L)) { stop("non-empty anchor in data that is not in margins") }
	tmatch <- all.indices[targets(data, id=TRUE)]
	if (any(tmatch==0L)) { stop("non-empty target in data that is not in margins") }
	
	return(data.frame(amatch, tmatch))
}	
