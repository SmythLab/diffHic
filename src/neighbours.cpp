#include "neighbours.h"

basic::basic(int w, int t, bool i, int x=0) : level(0), width(w), tlen(t), intra(i), exclude(x) {
	if (width < 0 || exclude < 0) { throw std::runtime_error("width values must be non-negative"); }
	if (exclude >= width) { throw std::runtime_error("exclusion width must be less than flank width"); }
}

void basic::restrain () {
	if (left < 0) { left=0; }
	if (intra) {
		if (right > row) { right=row+1; }
	} else if (right > tlen) { right=tlen; } // For intra's, right will hit diagonal; no need to worry about tlen.
	if (left > right) { left=right; } // Avoid negative areas in calling function.
}

bottomright::bottomright(int w, int t, bool i, int x) : basic(w, t, i, x) { level=-w; }

bool bottomright::bump_level () { 
	if (level >= 0) { return false; }
	++level;
	return true;
}

void bottomright::set(int a, int t) {
	row=a+level;
	left=(level < -exclude ? t : t+exclude+1);
	right=t+width+1; 
	restrain();
}

updown::updown(int w, int t, bool i, int x) : basic(w, t, i, x) { level=-w; }
	
bool updown::bump_level() { 
	if (level >= width) { return false; }
	++level;
	if (level==-exclude) { level=exclude+1; }
	return true; 
}

void updown::set(int a, int t) {
	row=a+level;
	left=t;
	right=t+1;
	restrain();
}

/* `leftright` is split into two regions, one for the left of the excluded zone
 * and another for the right of that zone.
 */

leftright1::leftright1(int w, int t, bool i, int x) : basic(w, t, i, x) {} 
	
bool leftright1::bump_level() { return false; }

void leftright1::set(int a, int t) { 
	row=a;
	left=t-width;
	right=t-exclude;
	restrain(); 
}	

leftright2::leftright2(int w, int t, bool i, int x) : basic(w, t, i, x) {} 

bool leftright2::bump_level() { return false; }

void leftright2::set(int a, int t) { 
	row=a;
	left=t+exclude+1;
	right=t+width+1;
	restrain(); 
}	

/* `allaround` is split into two rotationally-symmetric shapes; sort of
 * like the L-blocks in tetris, which - when fitted against each other
 * - form a 3x3 square with a hole in the middle. That hole is the 
 * excluded zone in this context, generalized for WxW squares.
 */

allaround1::allaround1(int w, int t, bool i, int x) : basic(w, t, i, x) { level=-w; } 

bool allaround1::bump_level () { 
	if (level >= width) { return false; }
	++level;
	return true;
}

void allaround1::set(int a, int t) {
	row=a+level;
	left=t-width;
	right=(level < -exclude ? t+exclude+1 : t-exclude);
	restrain();		
}

allaround2::allaround2(int w, int t, bool i, int x) : basic(w, t, i, x) { level=-w; }

bool allaround2::bump_level () { 
	if (level >= width) { return false; }
	++level;
	return true;
}

void allaround2::set(int a, int t) {
	row=a+level;
	left=(level > exclude ? t-exclude : t+exclude+1);
	right=t+width+1;
	restrain();		
}

