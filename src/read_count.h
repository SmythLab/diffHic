#ifndef READ_COUNT_H
#define READ_COUNT_H

#include "diffhic.h"

struct coord {
    coord (const int& a, const int& t, const int& l) : anchor(a), target(t), library(l) {}
	bool operator> (const coord& other) const {
		if (anchor > other.anchor) { return true; }
		else if (anchor == other.anchor) { return target>other.target; }
		else { return false; }
	}
	int anchor, target, library;
};

class binner {
public:
	binner(SEXP, SEXP, int, int);
	void fill(int*, bool*, std::deque<int>&);
	bool empty() const;
	int get_nlibs() const;
	int get_anchor() const;
private:
	const int fbin, lbin;
	const int* bptr;
	int nlibs;
	
	std::deque<const int*> aptrs, tptrs;
	std::deque<int> nums, indices;
	std::priority_queue<coord, std::deque<coord>, std::greater<coord> > next;
	
	int curab, curtb, curlib, curdex, lib;
	bool failed;
};

double nb_average(const int&, const int&, const double&, const double*, const int*, const double&);

#endif
