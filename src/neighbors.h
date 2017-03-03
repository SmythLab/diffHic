#ifndef NEIGHBORS_H
#define NEIGHBORS_H

#include "diffhic.h"

/* There are several virtual methods here:
 *
 * - set(), which specifies the bin pair of interest and computes the left/right boundaries of a given neighbourhood at a given anchor.
 *   For reference, the interaction space is defined here such that each row is an anchor value and each column is a target.
 *   (This is an upper-triangular space for intra-chromosomal interaction matrices.)
 *
 * - bump_level(), which changes the level of the anchor at which set() operates. 
 *   Possible levels are determined by the up/down boundaries of the neighbourhood.
 *
 * The idea is to have a set of bin pairs sorted by anchor and target indices.
 * At a given level, the algorithm runs through the bin pairs in order.
 * The abundance for the neighbourhood at this level is summed from all bin pairs at the same anchor and within the left/right boundaries.
 * This is repeated for all possible levels to get the total sum for the neighbourhood.
 *
 * Note that set() will happily give anchors that are negative or greater than the largest possible anchor.
 * These should be ignored, using 'continue' for negative anchors and 'break' for above-maximal anchors when looping across bin pairs.
 *
 */


struct basic {
	basic(int, int, bool, int);
	virtual void set (int, int)=0;
	virtual ~basic() {};
	virtual bool bump_level ()=0; 
	int row, left, right;
protected:
	int level, width, tlen;
	bool intra;
	int exclude;
	void restrain ();
};

struct bottomright : public basic { 
	bottomright(int, int, bool, int);
	~bottomright() {};
	bool bump_level();
	void set(int, int);
};

struct updown : public basic {
	updown(int, int, bool, int);
	~updown() {};
	bool bump_level();
	void set(int, int);
};

struct leftright : public basic {
	leftright(int, int, bool, int);
	~leftright() {};
	bool bump_level();
	void set(int, int);
protected:
    bool onleft;
};

enum Rangetype { FULLRANGE, LEFTSIDE, RIGHTSIDE };

struct allaround : public basic {
	allaround(int, int, bool, int);
	~allaround() {};
	bool bump_level();
	void set(int, int);
protected:
    Rangetype rangemode;
};

#endif
