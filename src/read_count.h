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

typedef std::priority_queue<coord, std::deque<coord>, std::greater<coord> > pair_queue;

void setup_pair_data (SEXP, const int, std::deque<const int*>&, std::deque<const int*>&, std::deque<int>&, std::deque<int>&);

class binner {
public:
	binner(SEXP, SEXP, int, int);
    ~binner();

	void fill();
	bool empty() const;
	int get_nlibs() const;
	int get_nbins() const;
	int get_anchor() const;

    const std::deque<int>& get_counts() const;
    const std::deque<int>& get_changed() const;
private:
	const int fbin, lbin, nbins;
	const int* bptr;
	int nlibs;

    std::deque<const int*> aptrs, tptrs;
	std::deque<int> nums, indices;
    pair_queue next;
	
	int curab, curtb, curlib, curdex, lib;
	bool failed;

    // Stuff that is visible to the calling class.
    std::deque<int> curcounts;
    std::deque<bool> ischanged;
    std::deque<int> waschanged;
};

#endif
