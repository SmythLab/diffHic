#include "diffhic.h"
//#define DEBUG 1

#ifdef DEBUG
#include <iostream>
#endif

/****************************************************************************
 ****************************************************************************
 * 
 * This function performs the cluster identification using a simple
 * single-linkage approach. Specifically, clusters are formed so long
 * as each component is no more than 'tol' away from at least one other
 * cluster component.
 *
 ****************************************************************************
 ****************************************************************************/

const int nothing=-1;

SEXP cluster_2d (SEXP start_a, SEXP start_t, SEXP end_a, SEXP end_t, SEXP tol, SEXP verbose) try {
	if (!isInteger(start_a) || !isInteger(start_t) || !isInteger(end_a) || !isInteger(end_t)) {
		throw std::runtime_error("anchor or target start and ends must be integer"); }
	const int npts=LENGTH(start_a);
	if (npts!=LENGTH(start_t) || npts!=LENGTH(end_a) || npts!=LENGTH(end_t)) { 
		throw std::runtime_error("lengths of coordinate vectors are not equal"); }
	if (!isInteger(tol) || LENGTH(tol)!=1) { 
		throw std::runtime_error("tolerance should be an integer scalar"); }
	if (!isLogical(verbose) || LENGTH(verbose)!=1) { 
		throw std::runtime_error("verbosity should be a logical scalar"); }

	const int width=asInteger(tol);
	const int* asptr=INTEGER(start_a),
		* aeptr=INTEGER(end_a),
		* tsptr=INTEGER(start_t),
		* teptr=INTEGER(end_t);
	const int verb=asLogical(verbose);
	
	// Setting up the output construct, which specifies the cluster ID of each point.
	SEXP output=PROTECT(allocVector(INTSXP, npts));
try {
	int* optr=INTEGER(output);
	
	/* We assume all points are sorted by start_a, start_t already. The idea is to
 	 * raster over the interaction space, assimilating all those which lie within 'tol'
 	 * of previous points. Overlaps with multiple points result in the formation of
 	 * synonyms which will be resolved at the end of the function.
 	 */
	std::map<int, std::pair<int, int> > landscape;
	std::map<int, std::pair<int, int> >::iterator itl, ito;
	std::set<std::pair<int, int> > synonyms;
	int ptdex=0, numids=0;
	bool beforestart;

	while (ptdex < npts) {
		const int& cura=asptr[ptdex];
		const int& curt=tsptr[ptdex];
		const int& endt=teptr[ptdex];
		const int& enda=aeptr[ptdex];
#ifdef DEBUG
		std::cout << "Adding: " << cura << " to " << enda <<", " << curt << " to " << endt << std::endl;
#endif		

		const int baset=curt-width;
		const int basea=cura-width;
		const int finisht=endt+width;
		int& myid=(optr[ptdex]=nothing);

		// Identifying the landscape element before or at the current target
		itl=landscape.lower_bound(baset);
		if (itl==landscape.end()) {
			; 
		} else {
			if ((itl->first)!=baset && itl!=landscape.begin()) {
				--itl;
				if ((itl->second).first > basea) { myid=(itl->second).second; }
			}
			
			// Running through everything in range, seeing if it overlaps anything else.
			while (itl!=landscape.end() && itl->first < finisht) {
 				if ((itl->second).first > basea) {
					const int& alternative=(itl->second).second;
					if (alternative==nothing) { 
						;
					} else if (myid==nothing) {
						myid=alternative;
					} else if (myid!=alternative) {
						if (alternative > myid) {
 	 					   	synonyms.insert(std::make_pair(myid, alternative));
						} else {
 	 					   	synonyms.insert(std::make_pair(alternative, myid));
						}
					}
 				}
				++itl;
			} 
		}
		if (myid==nothing) { // Incrementing to the next ID, if the current box doesn't overlap with anything.
			myid=numids;
			++numids;
		}

		// Adding the current endpoint of our box to the landscape.
		if (itl==landscape.begin()) {
 	 	   	itl=landscape.insert(itl, std::make_pair(endt, std::make_pair(nothing, nothing)));
		} else {
			--itl; 
			while (itl!=landscape.begin() && itl->first > endt) { --itl; } // Get to the first point in front of or equal to endt.

			/* We only need to change if it's different, remember; we want to keep the existing underlying element.
			 * We also only add it if the landscape is less than the added point, otherwise we'd be underground or redundant.
			 * Of course, the latter only applies if it's not a 'nothing' endpost.
			 */
			if (itl->first!=endt && ((itl->second).first < enda || (itl->second).second==nothing)) {
				const int& lastid=(itl->second).second;
				const int& lastheight=(itl->second).first;
				++itl;
				itl=landscape.insert(itl, std::make_pair(endt, std::make_pair(lastheight, lastid)));
			} else if (itl->first!=endt) {
				++itl; // Set to first element greater than endt, if nothing is equal to it. Ensures later `--itl' starts at the right point.
			}
		}

		/* After we find endt, curt is always smaller, so regardless of width, we can just traverse downwards.	
		 * We either delete the existing points if they are lower than the current height, or we alter it if 
		 * it was a 'nothing' post. Deletion is possible because there's no way that later elements with higher
		 * start anchors can overlap an end anchor lower than our current end anchor without also overlapping
		 * our current end anchor; moreover, if our current start anchor didn't overlap the deleted end anchor,
		 * then there's no way that later start anchors could overlap it either. So, the highest would be the best.
		 * Finally, we add a new start post if the preceding element is not overriding us.
		 */
		beforestart=(itl==landscape.begin());
		if (!beforestart) { 
			--itl;
			while (itl->first > curt) {
				if ((itl->second).second==nothing) {
					// Peeking to see if the previous element has a high 'enda'.
					--itl;
					if ((itl->second).first <= enda) {
						++itl;
						landscape.erase(itl--);
					} else {
						++itl;
						(itl->second).first=enda;
						(itl->second).second=myid;
						--itl; // landscape can't start with nothing, so don't worry about begin().
					}
				} else if ((itl->second).first <= enda) { 
					if (itl!=landscape.begin()) { 
						landscape.erase(itl--); 
					} else {
						landscape.erase(itl++);
						beforestart=true;
						break;
					}
				} else {
					if (itl==landscape.begin()) { 
						beforestart=true;
						break; 
					}
					--itl;
				}
 			}
		}
		if (beforestart) { // If the current target is before the landscape start, we must add it.
			itl=landscape.insert(itl, std::make_pair(curt, std::make_pair(enda, myid)));
		} else if ((itl->second).first < enda || (itl->second).second==nothing) { // Shouldn't be empty, after adding the endpoint, above. 
			if (itl->first==curt) {
				(itl->second).first=enda;
				(itl->second).second=myid;
			} else {
				++itl;
				itl=landscape.insert(itl, std::make_pair(curt, std::make_pair(enda, myid)));
			}
		}

		if (verb) { // Printing landscape.
			for (std::map<int, std::pair<int, int> >::const_iterator itx=landscape.begin(); itx!=landscape.end(); ++itx) {
                Rprintf("%i[%i](%i) ", itx->first, (itx->second).first, (itx->second).second);
			}
			Rprintf("\n");
		}

#ifdef DEBUG
		std::cout << "#### New landscape!" << std::endl;
		for (std::map<int, std::pair<int, int> >::const_iterator itx=landscape.begin(); itx!=landscape.end(); ++itx) {
			std::cout << itx->first << ": " << (itx->second).first << " (" << (itx->second).second << ")" << std::endl;
		}
#endif		
		++ptdex;
	}

	// Resolving synonyms using a recursive-ish algorithm. It's a bit slow but at least it'll work. I need to fill up the reverse, though.
	{ 
		std::set<std::pair<int, int> > alternate;
		for (std::set<std::pair<int, int> >::iterator its=synonyms.begin(); its!=synonyms.end(); ++its) { 
			alternate.insert(std::make_pair(its->second, its->first)); }
		synonyms.insert(alternate.begin(), alternate.end());
	}
#ifdef DEBUG
	std::cout << "Synonyms are: " << std::endl;
	for (std::set<std::pair<int, int> >::const_iterator itx=synonyms.begin(); itx!=synonyms.end(); ++itx) {
		std::cout << itx->first << "\t" << itx->second << std::endl;;
	}
	std::cout << "Iterating now!" << std::endl;
#endif
	
	std::deque<int> newids(numids, nothing);
	numids=1;
	while (!synonyms.empty()) {
		std::priority_queue<int, std::deque<int>, std::greater<int> > next;
		std::set<std::pair<int, int> >::iterator itx=synonyms.begin();
		int curid=itx->first;
		newids[curid]=numids;

		bool found=false;
		do {
			do {
				next.push(itx->second);
				newids[itx->second]=numids;
				synonyms.erase(itx++);
			} while (itx!=synonyms.end() && curid==itx->first);

			found=false;
			while (!next.empty()) {
				curid=next.top();
				itx=synonyms.lower_bound(std::make_pair(curid, nothing));
				do { next.pop(); } while (!next.empty() && next.top()==curid);

				if (itx!=synonyms.end() && itx->first==curid) {
					found=true;
					break;
				}
			}
		} while (found);
		++numids;
	}

	// Mopping up anything which doesn't have any synonyms.
	for (size_t i=0; i<newids.size(); ++i) {
		if (newids[i]==nothing) { 
			newids[i]=numids;
			++numids;
		}
	}
	
	for (int i=0; i<npts; ++i) { optr[i]=newids[optr[i]]; }
} catch (std::exception &e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

/****************************************************************************
 ****************************************************************************
 *
 * This function performs the cluster splitting, whereby the coordinates of 
 * each cluster set are computed, the subintervals identified, and then, 
 * everything is recomputed for each cluster subdivision. I've put this in
 * a separate function in order to let it breathe for a bit; otherwise, the 
 * clustering function itself becomes too large.
 *
 *
 ****************************************************************************
 ****************************************************************************/

SEXP split_clusters (SEXP id, SEXP start_a, SEXP start_t, SEXP end_a, SEXP end_t, SEXP maxw) try {
	if (!isInteger(id)) { throw std::runtime_error("cluster id values should be integer"); }
	const int npts=LENGTH(id);
	if (!isInteger(start_a) || !isInteger(start_t) || !isInteger(end_a) || !isInteger(end_t)) {
		throw std::runtime_error("anchor or target start and ends must be integer"); }
	if (npts!=LENGTH(start_a) || npts!=LENGTH(start_t) || npts!=LENGTH(end_a) || npts!=LENGTH(end_t)) { 
		throw std::runtime_error("lengths of coordinate vectors are not equal"); }
	if (!isInteger(maxw) || LENGTH(maxw)!=1) { 
		throw std::runtime_error("maximum width should be an integer scalar"); }

	const int width=asInteger(maxw);
	const int* asptr=INTEGER(start_a),
		* aeptr=INTEGER(end_a),
		* tsptr=INTEGER(start_t),
		* teptr=INTEGER(end_t),
		* iptr=INTEGER(id);

	// Getting the maximum ID, and constructing holding cells specifying the cluster starts and ends. 
	int maxid=0;	
	for (int i=0; i<npts; ++i) { 
		if (maxid < iptr[i]) { maxid=iptr[i]; }
	}
	++maxid;
	int* newas=(int*)R_alloc(maxid, sizeof(int));
	for (int i=0; i<maxid; ++i) { newas[i]=-1; }
	int* newts=(int*)R_alloc(maxid, sizeof(int));
	int* newae=(int*)R_alloc(maxid, sizeof(int));
	int* newte=(int*)R_alloc(maxid, sizeof(int));

	// Getting the extremes for the starts and ends.
	for (int i=0; i<npts; ++i) {
		const int& curid=iptr[i];
		if (newas[curid]==-1) { 
			newas[curid]=asptr[i];
			newts[curid]=tsptr[i];
			newae[curid]=aeptr[i];
			newte[curid]=teptr[i];
		} else {
			if (asptr[i] < newas[curid]) { newas[curid]=asptr[i]; }
			if (tsptr[i] < newts[curid]) { newts[curid]=tsptr[i]; }
			if (aeptr[i] > newae[curid]) { newae[curid]=aeptr[i]; }
			if (teptr[i] > newte[curid]) { newte[curid]=teptr[i]; }
		}
	}

	// Going through each cluster and filling in the subinterval width.
	double* newsuba=(double*)R_alloc(maxid, sizeof(double));
	double* newsubt=(double*)R_alloc(maxid, sizeof(double));
	int * newnumt =(int*)R_alloc(maxid, sizeof(int));
	int * newidstart=(int*)R_alloc(maxid, sizeof(int));
	double diff;
	int mult, prod, last=1;

	for (int i=0; i<maxid; ++i) {
		if (newas[i]==-1) { continue; }
		diff=newae[i]-newas[i];
		if (diff > width) { 
			mult=int(diff/width+0.5);
			prod=mult;
			newsuba[i]=diff/mult;
		} else { 
			newsuba[i]=-1; 
			prod=1;
		}

		diff=newte[i]-newts[i];
		if (diff > width) { 
			newnumt[i]=int(diff/width+0.5);
			prod*=newnumt[i];
			newsubt[i]=diff/newnumt[i];
		} else {
 			newnumt[i]=1;
			newsubt[i]=-1; 
		}

		newidstart[i]=last;
		last+=prod;
	}

	/* Now, going through the original points and figuring out which subparition it fits into.
 	 * We use the midpoint of the original bins to decide if something should go in one place
 	 * or the other.
 	 */
	SEXP output=PROTECT(allocVector(INTSXP, npts));
	try {
		int* optr=INTEGER(output);
		for (int i=0; i<npts; ++i) {
			const int& curid=iptr[i];
			if (newas[curid]==-1) { continue; }
			int& extra=(optr[i]=newidstart[curid]);
			if (newsuba[curid]>0) {
				extra+=int( (0.5*double(aeptr[i]+asptr[i])-newas[curid])/newsuba[curid] ) * newnumt[curid];
			}
			if (newsubt[curid]>0) {
				extra+=int( (0.5*double(teptr[i]+tsptr[i])-newts[curid])/newsubt[curid] );
			} 
		}
	} catch (std::exception& e){
		UNPROTECT(1);
		throw;
	}

	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

/*************************************************************
 * This function performs the calculation of the bounding box
 * for each cluster.
 *************************************************************/

SEXP get_bounding_box (SEXP ids, SEXP starts, SEXP ends) try {
	if (!isInteger(ids)) { throw std::runtime_error("ID vector should be integer"); }
	if (!isInteger(starts)) { throw std::runtime_error("start vector should be integer"); }
	if (!isInteger(ends)) { throw std::runtime_error("end vector should be integer"); }
	
	const int npts=LENGTH(ids);
	if (LENGTH(starts)!=npts || LENGTH(ends)!=npts) { throw std::runtime_error("vectors are not of same length"); }
	const int * iptr=INTEGER(ids);
	const int * sptr=INTEGER(starts);
	const int * eptr=INTEGER(ends);

	int maxid=0;
	for (int i=0; i<npts; ++i) { 
		if (iptr[i] > maxid) { maxid=iptr[i]; }
	}

	SEXP output=PROTECT(allocVector(VECSXP, 3));
try {
	SET_VECTOR_ELT(output, 0, allocVector(INTSXP, maxid));	
	int* first_ptr=INTEGER(VECTOR_ELT(output, 0));
	for (int i=0; i<maxid; ++i) { first_ptr[i] = -1; }
	SET_VECTOR_ELT(output, 1, allocVector(INTSXP, maxid));	
	int* start_ptr=INTEGER(VECTOR_ELT(output, 1));
	SET_VECTOR_ELT(output, 2, allocVector(INTSXP, maxid));	
	int* end_ptr=INTEGER(VECTOR_ELT(output, 2));

	// To deal with 1-based indexing.
	--start_ptr;
	--end_ptr;
	--first_ptr;
	for (int i=0; i<npts; ++i) { 
		const int& current=iptr[i];
		if (first_ptr[current]==-1) {
			first_ptr[current]=i+1;
			start_ptr[current]=sptr[i];
			end_ptr[current]=eptr[i];
		} else { 
			if (start_ptr[current] > sptr[i]) { start_ptr[current]=sptr[i]; } 
			if (end_ptr[current] < eptr[i]) { end_ptr[current]=eptr[i]; }
		}
	}	

	++first_ptr;
	for (int i=0; i<maxid; ++i) { 
		if (first_ptr[i]==-1) { throw std::runtime_error("missing entries in the ID vector"); }
	}
} catch (std::exception& e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception &e) {
	return mkString(e.what());
}

