#include "read_count.h"
#include "neighbours.h"
#include <cmath>

const double low_value=std::pow(10.0, -10.0);
const double MULT=1000000;
typedef std::pair<int, int> floater;

floater nb_average(const int& nlibs, const int& maxit, const double& tolerance, const double* offset, 
		const int* y, const double& disp) { // Code copied from glm_one_group in edgeR's src.

	bool nonzero=false;
	double cur_beta=0;
	for (int j=0; j<nlibs; ++j) {
		const int& cur_val=y[j];
		if (cur_val>low_value) {
			cur_beta+=cur_val/std::exp(offset[j]);
			nonzero=true;
		}
	}
	cur_beta=std::log(cur_beta/nlibs);
	if (!nonzero) { return std::make_pair(0, 0); }

	double dl, info, mu, step, denominator;
	for (int i=0; i<maxit; ++i) {
		dl=0;
 	    info=0;
		for (int j=0; j<nlibs; ++j) {
			mu=std::exp(cur_beta+offset[j]), denominator=1+mu*disp;
			dl+=(y[j]-mu)/denominator;
			info+=mu/denominator;
		}
		step=dl/info;
		cur_beta+=step;
		if (std::abs(step)<tolerance) { break; }
	}

	// Getting integer values.
	cur_beta=std::exp(cur_beta);
	int int_comp=int(cur_beta);
	int dec_comp=int((cur_beta-double(int_comp))*MULT);
	return std::make_pair(int_comp, dec_comp);
}

/*********************************************************************
 ******************** Main looping. **********************************
 *********************************************************************/


SEXP count_background(SEXP all, SEXP bin, SEXP back_width, SEXP filter, 
		SEXP first_target_bin, SEXP last_target_bin, SEXP first_anchor_bin, SEXP last_anchor_bin,
		SEXP max_it, SEXP tolerance, SEXP offsets, SEXP dispersion,
		SEXP prior_count) try {
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int f=asInteger(filter);
	if (!isInteger(back_width) || LENGTH(back_width)!=1) { throw std::runtime_error("width of neighbourhood regions must be an integer scalar"); }
	const int bwidth=asInteger(back_width);

	// Getting the indices of the first and last bin on the target and anchor chromosomes.
	if (!isInteger(first_target_bin) || LENGTH(first_target_bin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int ftbin=asInteger(first_target_bin);
	if (!isInteger(last_target_bin) || LENGTH(last_target_bin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int ltbin=asInteger(last_target_bin);
	if (!isInteger(first_anchor_bin) || LENGTH(first_anchor_bin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int fabin=asInteger(first_anchor_bin);
	if (!isInteger(last_anchor_bin) || LENGTH(last_anchor_bin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int labin=asInteger(last_anchor_bin);
	const bool intra=(fabin==ftbin);

	// Setting up the binning engine.
	binner engine(all, bin, ftbin, ltbin);
	const int nlibs=engine.get_nlibs();

	// Setting up the NB average stuff.
	if (!isInteger(max_it) || LENGTH(max_it)!=1) { throw std::runtime_error("maximum number of iterations must be an integer scalar"); }
	const int maxit=asInteger(max_it);
	if (!isReal(tolerance) || LENGTH(tolerance)!=1) { throw std::runtime_error("tolerance must be a double-precision vector"); }
	const double tol=asReal(tolerance);
	if (!isReal(offsets) || LENGTH(offsets)!=nlibs) { throw std::runtime_error("offsets must be a double-precision vector of length equal to number of libraries"); }
	const double* offptr=REAL(offsets);
	if (!isReal(dispersion) || LENGTH(dispersion)!=1) { throw std::runtime_error("dispersion must be a double-precision vector"); }
	const double disp=asReal(dispersion);
	if (!isReal(prior_count) || LENGTH(prior_count)!=1) { throw std::runtime_error("prior count must be a double-precision vector"); }
	const double prior=asReal(prior_count);

	// Setting up the memory containers.
	const int ntbins=ltbin-ftbin+1, nabins=labin-fabin+1;
	int* curcounts=(int*)R_alloc(nlibs*ntbins, sizeof(int)); 
	bool* ischanged=(bool*)R_alloc(ntbins, sizeof(bool));
	for (int i=0; i<ntbins; ++i) { ischanged[i]=false; }
	std::deque<int> waschanged, counts, anchors, targets;
	std::deque<floater> averages;

	// Stuff required to compute the neighbourhood.
	std::deque<int> ref_targets, ref_anchors;
	std::deque<floater> ref_ave;
	size_t nmodes=(intra ? 4 : 3), mode;
	std::deque<std::deque<floater> > neighbourave(nmodes);
    std::deque<std::deque<int> > neighbourarea(nmodes);
	basic * base_ptr=NULL;
	size_t saved_dex=0; // Points to the entry in 'anchors' to be interrogated.
	size_t n_to_drop;

	// Other assorted sundries.
	size_t vecdex;
	int rowdex, countsum, curlib, curanchor;
	floater current_average;
	int leftbound, rightbound, leftdex, rightdex, desired_anchor;
   	size_t saved_copy_dex;

	while (1) { 
		if (!engine.empty()) {
			engine.fill(curcounts, ischanged, waschanged);
			curanchor=engine.get_anchor() - fabin;

			// Adding it to the main list, if it's large enough.
			for (vecdex=0; vecdex<waschanged.size(); ++vecdex) {
				ref_anchors.push_back(curanchor);
				ref_targets.push_back(waschanged[vecdex]);

				rowdex=waschanged[vecdex]*nlibs;
				current_average=nb_average(nlibs, maxit, tol, offptr, curcounts+rowdex, disp);
				ref_ave.push_back(current_average);

				countsum=0;
				for (curlib=0; curlib<nlibs; ++curlib) { countsum+=curcounts[rowdex+curlib]; }
				if (countsum >= f) { 
					anchors.push_back(curanchor);
					targets.push_back(ref_targets.back());
					averages.push_back(current_average);
					for (curlib=0; curlib<nlibs; ++curlib) { counts.push_back(curcounts[rowdex+curlib]); }
				}
				// Resetting the ischanged vector for the next stretch of bin anchors.
				ischanged[waschanged[vecdex]]=false;
			}
			waschanged.clear();

			// Matching size.
			for (mode=0; mode<nmodes; ++mode) { 
				neighbourarea[mode].resize(anchors.size());
				neighbourave[mode].resize(averages.size()); 
			} 
		}

		if (saved_dex < anchors.size() && (engine.empty() || ref_anchors.back() >= anchors[saved_dex] + bwidth)) { 
			bottomright br(bwidth, ntbins, intra);
			leftright lr(bwidth, ntbins, intra);
			updown ud(bwidth, ntbins, intra);
			allaround aa(bwidth, ntbins, intra);
			const int& saved_anchor=anchors[saved_dex];
			int fullsize=ref_anchors.size();

			// Computing the neighbourhood count for all bin pairs with anchors of 'anchors[saved_dex]'.
			for (mode=0; mode<nmodes; ++mode) { 
				switch (mode) {
					case 0: base_ptr=&lr; break;
					case 1: base_ptr=&ud; break;
					case 2: base_ptr=&aa; break;
					case 3: base_ptr=&br; break;
				}

				leftdex=rightdex=0;
				current_average.first=current_average.second=0;
				do { 
					for (saved_copy_dex=saved_dex; saved_copy_dex < anchors.size() && 
							anchors[saved_copy_dex]==saved_anchor; ++saved_copy_dex) { 
						base_ptr->set(saved_anchor, targets[saved_copy_dex]);
						desired_anchor=base_ptr->row;
						if (desired_anchor < 0 || desired_anchor >= nabins) { break; }
						leftbound=base_ptr->left;
						rightbound=base_ptr->right;

						// Identifying all reference elements associated with the neighbourhood of this saved_copy_dex at the current level.
						while (rightdex < fullsize && (ref_anchors[rightdex] < desired_anchor || 
									(ref_anchors[rightdex]==desired_anchor && ref_targets[rightdex] < rightbound))) { 
							current_average.first+=ref_ave[rightdex].first;
							current_average.second+=ref_ave[rightdex].second;
							++rightdex;
						}
						while (leftdex < fullsize && (ref_anchors[leftdex] < desired_anchor || 
									(ref_anchors[leftdex]==desired_anchor && ref_targets[leftdex] < leftbound))) {
							current_average.first-=ref_ave[leftdex].first;
							current_average.second-=ref_ave[leftdex].second;
							++leftdex;
						}
		
						neighbourave[mode][saved_copy_dex].first+=current_average.first;
						neighbourave[mode][saved_copy_dex].second+=current_average.second;
						neighbourarea[mode][saved_copy_dex]+=rightbound-leftbound;
						if (base_ptr->discard_self()) { 
							neighbourave[mode][saved_copy_dex].first-=averages[saved_copy_dex].first;
							neighbourave[mode][saved_copy_dex].second-=averages[saved_copy_dex].second;
							--(neighbourarea[mode][saved_copy_dex]);
						}
					}
					/* Each bump_level will result in an increase in desired_anchor, so we can hot-start 
					 * from the previous left/rightdex. Both indices MUST increase (even when saved_copy_dex 
					 * is reset to saved_dex) as they've been stuck on the previous desired_anchor.
					 */
				} while (base_ptr->bump_level());
			}
			while (saved_dex < anchors.size() && anchors[saved_dex]==saved_anchor) { ++saved_dex; } // Shifting onwards.
		}

		// Dropping elements in the reference data that will no longer be used, to save memory.
		if (saved_dex < anchors.size()) { 
			n_to_drop=0;
			while (n_to_drop < ref_anchors.size() && ref_anchors[n_to_drop] < anchors[saved_dex] - bwidth) { ++n_to_drop; }
			ref_anchors.erase(ref_anchors.begin(), ref_anchors.begin()+n_to_drop);
			ref_targets.erase(ref_targets.begin(), ref_targets.begin()+n_to_drop);
			ref_ave.erase(ref_ave.begin(), ref_ave.begin()+n_to_drop);
		} else if (engine.empty()) { break; }
	}

	// Storing into R objects.
	SEXP output=PROTECT(allocVector(VECSXP, 4));
	try {
		const int ncombos=anchors.size();
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, ncombos));
		int* aoptr=INTEGER(VECTOR_ELT(output, 0));
		SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ncombos));
		int* toptr=INTEGER(VECTOR_ELT(output, 1));
		SET_VECTOR_ELT(output, 2, allocMatrix(INTSXP, ncombos, nlibs));
		SET_VECTOR_ELT(output, 3, allocVector(REALSXP, ncombos));
		double* foptr=REAL(VECTOR_ELT(output, 3));
	
		std::deque<int*> coptrs(nlibs);
		for (curlib=0; curlib<nlibs; ++curlib) {
			if (curlib==0) { coptrs[curlib]=INTEGER(VECTOR_ELT(output, 2)); } 
			else { coptrs[curlib]=coptrs[curlib-1]+ncombos; }
		}	
			
		// Iterating across and filling both the matrix and the components.
		int cdex=-1;
		double backaverage, maxaverage;
		for (vecdex=0; vecdex<anchors.size(); ++vecdex) {
			aoptr[vecdex]=anchors[vecdex]+fabin;
			toptr[vecdex]=targets[vecdex]+ftbin;
			for (curlib=0; curlib<nlibs; ++curlib) {
				++cdex; 			   	
				coptrs[curlib][vecdex]=counts[cdex]; 
			}
			
			foptr[vecdex]=double(averages[vecdex].first)+double(averages[vecdex].second)/MULT;
			maxaverage=0;
			for (mode=0; mode<nmodes; ++mode) {
				backaverage=double(neighbourave[mode][vecdex].first)+double(neighbourave[mode][vecdex].second)/MULT;
				if (neighbourarea[mode][vecdex]) { 
					backaverage/=neighbourarea[mode][vecdex];
					if (backaverage > maxaverage) { maxaverage=backaverage; }
				}
			}
			foptr[vecdex]+=prior;
			foptr[vecdex]/=maxaverage+prior;			
		}
	} catch (std::exception& e) { 
		UNPROTECT(1);
		throw;
	}

	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}


