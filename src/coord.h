#ifndef COORD_H
#define COORD_H

struct coord {
    coord (const int& a, const int& t, const int& l) : anchor(a), target(t), library(l) {}
	bool operator> (const coord& other) const {
		if (anchor > other.anchor) { return true; }
		else if (anchor == other.anchor) { return target>other.target; }
		else { return false; }
	}
	int anchor, target, library;
};

#endif

