#include "read_count.h"
#include <cmath>

const double low_value=std::pow(10.0, -10.0);

double nb_average(const int& nlibs, const int& maxit, const double& tolerance, const double* offset, 
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
	if (!nonzero) { return 0; }

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
	return std::exp(cur_beta);
}

