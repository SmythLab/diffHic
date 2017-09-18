#include "diffhic.h"
#include "utils.h"
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

SEXP cluster_2d (SEXP start_anchor1, SEXP start_anchor2, SEXP end_anchor1, SEXP end_anchor2, SEXP tol, SEXP verbose) {
    BEGIN_RCPP
    Rcpp::IntegerVector as1(start_anchor1), as2(start_anchor2), ae1(end_anchor1), ae2(end_anchor2);
    const int npts=as1.size();
    if (npts!=as2.size() || npts!=ae1.size() || npts!=ae2.size()) { 
        throw std::runtime_error("lengths of coordinate vectors are not equal"); 
    }

    const int width=check_integer_scalar(tol, "tolerance");
    const bool verb=check_logical_scalar(verbose, "verbosity");

	// Setting up the output construct, which specifies the cluster ID of each point.
    Rcpp::IntegerVector output(npts, nothing);
	
	/* We assume all points are sorted by start_a1, start_a2 already. The idea is to
 	 * raster over the interaction space, assimilating all those which lie within 'tol'
 	 * of previous points. Overlaps with multiple points result in the formation of
 	 * synonyms which will be resolved at the end of the function.
 	 */
	std::map<int, std::pair<int, int> > landscape;
	std::set<std::pair<int, int> > synonyms;
	int ptdex=0, numids=0;

    while (ptdex < npts) {
        const int& cura=as1[ptdex];
        const int& enda=ae1[ptdex];
        const int& curt=as2[ptdex];
        const int& endt=ae2[ptdex];
#ifdef DEBUG
		std::cout << "Adding: " << cura << " to " << enda <<", " << curt << " to " << endt << std::endl;
#endif		

		const int baset=curt-width;
		const int basea=cura-width;
		const int finisht=endt+width;
		int& myid=output[ptdex];

		// Identifying the landscape element before or at the current target
		auto itl=landscape.lower_bound(baset);
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
		bool beforestart=(itl==landscape.begin());
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
			for (auto itx=landscape.begin(); itx!=landscape.end(); ++itx) {
                Rprintf("%i[%i](%i) ", itx->first, (itx->second).first, (itx->second).second);
			}
			Rprintf("\n");
		}

#ifdef DEBUG
		std::cout << "#### New landscape!" << std::endl;
		for (auto itx=landscape.begin(); itx!=landscape.end(); ++itx) {
			std::cout << itx->first << ": " << (itx->second).first << " (" << (itx->second).second << ")" << std::endl;
		}
#endif		
		++ptdex;
	}

	// Resolving synonyms using a recursive-ish algorithm, which is a bit slow but at least it'll work. 
    // I need to fill up the reverse links, though, to guarantee that it will work properly.
	{ 
		std::set<std::pair<int, int> > alternate;
		for (const auto& link : synonyms) { 
			alternate.insert(std::make_pair(link.second, link.first)); 
        }
		synonyms.insert(alternate.begin(), alternate.end());
	}
#ifdef DEBUG
	std::cout << "Synonyms are: " << std::endl;
	for (const auto& link : synonyms) { 
		std::cout << link.first << "\t" << link.second << std::endl;;
	}
	std::cout << "Iterating now!" << std::endl;
#endif
	
	std::vector<int> newids(numids, nothing);
	numids=1;
	while (!synonyms.empty()) {
		std::priority_queue<int, std::deque<int>, std::greater<int> > next;
		auto itx=synonyms.begin();
		int curid=itx->first;
		newids[curid]=numids;

		bool found=false;
		do {
			do {
                const int& alt=itx->second;
				next.push(alt);
				newids[alt]=numids;
				synonyms.erase(itx++);
			} while (itx!=synonyms.end() && curid==itx->first);

			found=false;
			while (!next.empty()) {
				curid=next.top();
				do { 
                    next.pop(); 
                } while (!next.empty() && next.top()==curid);

				itx=synonyms.lower_bound(std::make_pair(curid, nothing));
				if (itx!=synonyms.end() && itx->first==curid) {
					found=true;
					break;
				}
			}
		} while (found);
		++numids;
	}

	// Mopping up anything which doesn't have any synonyms.
	for (auto& id : newids) { 
		if (id==nothing) { 
			id=numids;
			++numids;
		}
	}
	
	for (auto& o : output) { 
        o=newids[o];
    }
    return output;
    END_RCPP
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

SEXP split_clusters(SEXP clustid, SEXP start_anchor1, SEXP start_anchor2, SEXP end_anchor1, SEXP end_anchor2, SEXP maxw) {
    BEGIN_RCPP
    Rcpp::IntegerVector id(clustid), as1(start_anchor1), as2(start_anchor2), ae1(end_anchor1), ae2(end_anchor2);
    const int npts=id.size();
    if (npts!=as1.size() || npts!=as2.size() || npts!=ae1.size() || npts!=ae2.size()) { 
        throw std::runtime_error("lengths of coordinate vectors are not equal"); 
    }
    const int width=check_integer_scalar(maxw, "maximum width");
   
	// Getting the maximum ID, and constructing holding cells specifying the cluster starts and ends. 
    const int maxid=(npts ? (*std::max_element(id.begin(), id.end())) + 1 : 0);
    std::vector<int> newas1(maxid, nothing), newas2(maxid), newae1(maxid), newae2(maxid);

	// Getting the extremes for the starts and ends.
	for (int i=0; i<npts; ++i) {
		const int& curid=id[i];
		if (newas1[curid]==nothing) { 
			newas1[curid]=as1[i];
			newas2[curid]=as2[i];
			newae1[curid]=ae1[i];
			newae2[curid]=ae2[i];
		} else {
			if (as1[i] < newas1[curid]) { newas1[curid]=as1[i]; }
			if (as2[i] < newas2[curid]) { newas2[curid]=as2[i]; }
			if (ae1[i] > newae1[curid]) { newae1[curid]=ae1[i]; }
			if (ae2[i] > newae2[curid]) { newae2[curid]=ae2[i]; }
		}
	}

	// Going through each cluster and filling in the subinterval width.
    std::vector<double> newsuba1(maxid, -1), newsuba2(maxid, -1);
    std::vector<int> newnum(maxid, 1), newidstart(maxid);
	int last=1;

	for (int i=0; i<maxid; ++i) {
		if (newas1[i]==nothing) { 
            continue; 
        }
		double diff=newae1[i]-newas1[i];
        int prod=1;
		if (diff > width) { 
			const int mult=int(diff/width+0.5);
			prod=mult;
			newsuba1[i]=diff/mult;
		}

		diff=newae2[i]-newas2[i];
		if (diff > width) { 
            const int mult=int(diff/width+0.5);
			newnum[i]=mult;
			prod*=mult;
			newsuba2[i]=diff/mult;
		}

		newidstart[i]=last;
		last+=prod;
	}

	/* Now, going through the original points and figuring out which subparition it fits into.
 	 * We use the midpoint of the original bins to decide if something should go in one place
 	 * or the other.
 	 */
    Rcpp::IntegerVector output(npts);
    for (int i=0; i<npts; ++i) {
        const int& curid=id[i];
		int& extra=(output[i]=newidstart[curid]);
        if (newsuba1[curid]>0) {
			extra+=int( (0.5*double(ae1[i]+as1[i])-newas1[curid])/newsuba1[curid] ) * newnum[curid];
		}
		if (newsuba2[curid]>0) {
			extra+=int( (0.5*double(ae2[i]+as2[i])-newas2[curid])/newsuba2[curid] );
        } 
    }
	return output;
    END_RCPP
}

/*************************************************************
 * This function performs the calculation of the bounding box
 * for each cluster.
 *************************************************************/

SEXP get_bounding_box (SEXP ids, SEXP starts, SEXP ends) {
    BEGIN_RCPP
    Rcpp::IntegerVector _ids(ids), _starts(starts), _ends(ends);	
	const int npts=_ids.size();
	if (_starts.size()!=npts || _ends.size()!=npts) { 
        throw std::runtime_error("vectors are not of same length"); 
    }

	const int maxid=(npts ? *std::max_element(_ids.begin(), _ids.end()) : 0);
    Rcpp::IntegerVector out_first(maxid, nothing), out_start(maxid), out_end(maxid);

	for (int i=0; i<npts; ++i) { 
		const int current=_ids[i]-1; // 1-based indexing.
		if (out_first[current]==nothing) {
			out_first[current]=i+1;
			out_start[current]=_starts[i];
			out_end[current]=_ends[i];
		} else { 
			if (out_start[current] > _starts[i]) { out_start[current]=_starts[i]; } 
			if (out_end[current] < _ends[i]) { out_end[current]=_ends[i]; }
		}
	}	

    for (const auto& f : out_first) { 
		if (f==nothing) { throw std::runtime_error("missing entries in the ID vector"); }
	}
	return Rcpp::List::create(out_first, out_start, out_end);
    END_RCPP
}

