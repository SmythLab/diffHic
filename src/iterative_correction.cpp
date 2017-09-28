#include "diffhic.h"
#include "utils.h"

SEXP iterative_correction(SEXP avecount, SEXP anchor1, SEXP anchor2, SEXP local, 
		SEXP nfrags, SEXP iter, SEXP exlocal, SEXP lowdiscard, SEXP winsorhigh) {
    BEGIN_RCPP

	// Checking vector type and length.
    const Rcpp::NumericVector ac(avecount);
    const Rcpp::IntegerVector a1(anchor1), a2(anchor2);
    const Rcpp::LogicalVector l(local);
    const int npairs=ac.size();
    if (npairs!=a1.size() || npairs!=a2.size() || npairs!=l.size()) { 
        throw std::runtime_error("lengths of input vectors are not equal"); 
    }

	// Checking scalar type.
    const int numfrags=check_integer_scalar(nfrags, "number of fragments");
    const int iterations=check_integer_scalar(iter, "number of iterations");
    const int excluded=check_integer_scalar(exlocal, "exclusion specifier");
	const double discarded=check_numeric_scalar(lowdiscard, "proportion to discard");
    const double winsorized=check_numeric_scalar(winsorhigh, "proportion to winsorize");
	
	// Setting up a SEXP to hold the working probabilities (ultimately the truth).
    Rcpp::NumericVector working(ac.begin(), ac.end());

	// Excluding highly local interactions.
	int num_values=npairs;
	bool drop_intras=(excluded==NA_INTEGER);
	if (drop_intras || excluded >= 0) {
		for (int j=0; j<npairs; ++j) {
			if (l[j] && (drop_intras || a1[j]-a2[j] <= excluded)) { 
				working[j]=R_NaReal; 
				--num_values;
			}
		}
	}

	/**************************************************
 	 * Getting rid of unstable or extreme elements.
 	 **************************************************/

	// Winsorizing for the most abundant interactions.
	const int todrop=int(double(num_values)*winsorized);
	if (todrop>0) {
        std::vector<std::pair<double, int> > ordered(npairs);
        for (int i=0; i<npairs; ++i) {
            ordered[i].first=ac[i];
            ordered[i].second=i;
        }
        std::sort(ordered.begin(), ordered.end());
	
        // Identifying the winsorizing value.		
        int curcount=0;
        double winsor_val=R_NaReal;
        for (int i=npairs-1; i>=0; --i) { 
            const int& curo=ordered[i].second;
            const double& val=working[curo];
            if (ISNA(val)) { continue; }
            ++curcount;
            if (curcount > todrop) { 
                winsor_val=val;
                break; 
            }
        }

		// Applying said value to all high-abundance interactions.
	    if (ISNA(winsor_val)) { throw std::runtime_error("specified winsorizing proportion censors all data"); }
		curcount=0;
        for (int i=npairs-1; i>=0; --i) { 
            const int& curo=ordered[i].second;
            double& val=working[curo];
            if (ISNA(val)) { continue; }
            val=winsor_val;
            ++curcount;
            if (curcount > todrop) { break; }
        }		
	}

	// Bias and coverage vectors.
    Rcpp::NumericVector out_bias(numfrags, R_NaReal);
	for (int pr=0; pr<npairs; ++pr) { 
		if (!ISNA(working[pr])) { 
            out_bias[a1[pr]-1]=1;
            out_bias[a2[pr]-1]=1; 
        }
	}
    std::vector<double> coverage(numfrags);

	// Removing low-abundance fragments.
	if (discarded > 0) {
		for (int pr=0; pr<npairs; ++pr) {  // Computing the coverage.
            const double& wprob=working[pr];
			if (ISNA(wprob)) { continue; }
			coverage[a1[pr]-1]+=wprob;
			coverage[a2[pr]-1]+=wprob;
        }

        std::vector<std::pair<double, int> > ordering(numfrags);
		for (int i=0; i<numfrags; ++i) { 
            ordering[i].first=coverage[i];
            ordering[i].second=i; 
        }
		std::sort(ordering.begin(), ordering.end());
		
        // Ignoring empty rows.
	    int counter=0;
	    while (counter < numfrags && ISNA(out_bias[ordering[counter].second])){ ++counter; }
				
		// Filtering out low-abundance rows.
	    const int leftover=int(discarded*double(numfrags-counter))+counter;
		while (counter < leftover) { 
            out_bias[ordering[counter].second]=R_NaReal;
		    ++counter; 
		}

		// Setting everything else back to zero.
        while (counter < numfrags) { 
            coverage[ordering[counter].second]=0;
            ++counter;
        }
	
		// Propogating the filters through the interactions to remove anything with an anchor in the removed fragments.
		for (int pr=0; pr<npairs; ++pr) {
			if (ISNA(working[pr])) { continue; }
			if (ISNA(out_bias[a1[pr]-1]) || ISNA(out_bias[a2[pr]-1])) { 
                working[pr]=R_NaReal; 
            }  
		}

        // Resetting the values.
        std::fill(out_bias.begin(), out_bias.end(), R_NaReal);
        for (int pr=0; pr<npairs; ++pr) { 
            if (!ISNA(working[pr])) { 
                out_bias[a1[pr]-1]=1;
                out_bias[a2[pr]-1]=1; 
            }
        }
	}
	
    /********************************************************
 	 * Now, actually performing the iterative correction. ***
 	 ********************************************************/

    // Something to hold a diagnostic (i.e., the maximum step).
    Rcpp::NumericVector diagnostic(iterations);

	for (int it=0; it<iterations; ++it) {
		for (int pr=0; pr<npairs; ++pr) {  // Computing the coverage (ignoring locals, if necessary).
            const double& wprob=working[pr];
			if (ISNA(wprob)) { continue; }
			coverage[a1[pr]-1]+=wprob;
			coverage[a2[pr]-1]+=wprob;
        }
		
		/* Computing the 'aditional bias' vector, to take the mean across all fragments.
		 * The exact value of the mean doesn't really matter, it just serves to slow down 
		 * the algorithm to avoid overshooting the mark and lack of convergence. Everything
		 * ends up being the same so long as the column sums approach 1.
		 */
//		double mean_of_all=0;
//		for (int i=0; i<numfrags; ++i) { if (!ISNA(coverage[i])) { mean_of_all+=coverage[i]; } }
//		mean_of_all/=numfrags;
//		for (int i=0; i<numfrags; ++i) { if (!ISNA(coverage[i])) { coverage[i]/=mean_of_all; } }

		/* Okay, that didn't really work. The strategy that is described in the paper fails to
 		 * give me contact probabilities that sum to 1. So instead, I'm using the coverage
 		 * itself to divide the working probabilities. Square rooting is necessary to avoid
 		 * instability during iteration.
 		 */
		for (int i=0; i<numfrags; ++i) { 
            if (!ISNA(out_bias[i])) { 
                coverage[i]=std::sqrt(coverage[i]); 
            } 
        }

		// Dividing the working matrix with the (geometric mean of the) additional biases.	
		for (int pr=0; pr<npairs; ++pr){ 
            double& wprob=working[pr];
			if (!ISNA(wprob)) { 
                wprob/=coverage[a1[pr]-1] * coverage[a2[pr]-1]; 
            }
        }
		
		// Multiplying the biases by additional biases, and cleaning out coverage. We store
		// the maximum step to see how far off convergence it is.
		double& maxed=diagnostic[it];
		for (int i=0; i<numfrags; ++i) {			
            double& ob=out_bias[i];
			if (!ISNA(ob)) {
				double& cur_cov=coverage[i];
				if (cur_cov > maxed) { 
                    maxed=cur_cov; 
                } else if (1/cur_cov > maxed) { 
                    maxed=1/cur_cov; 
                }
                ob*=cur_cov;
				cur_cov=0;
			}
		}
    }

	/* Recalculating the contact probabilities, using the estimated biases, 
	 * to get normalized values for Winsorized or discarded bin pairs.
	 * Discarded bins can't be salvaged, though.
	 */
	for (int pr=0; pr<npairs; ++pr) {  
        const double& ob1=out_bias[a1[pr]-1];
        const double& ob2=out_bias[a2[pr]-1];
		if (ISNA(ob1) || ISNA(ob2)) { continue; }
		working[pr] = ac[pr]/(ob1*ob2);
	}
    
	return Rcpp::List::create(working, out_bias, diagnostic);
    END_RCPP
}

