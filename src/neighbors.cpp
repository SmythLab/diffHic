#include "neighbors.h"

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

/* The bottomright class identifies all bin pairs in a square with sides 'w'.
 * The bin pair of interest lies at the topleft corner.
 */

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

/* The updown class identifies all bin pairs in a vertical line of length 'w*2+1'.
 * The bin pair of interest lies in the centre.
 */

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

/* The leftright class identifies all bin pairs in a horizontal line of length 'w*2+1'.
 * The bin pair of interest lies in the centre. Here, bumping is done to cycle twice over 
 * each level; once to get the left side of the line, and again to get the right side.
 */

leftright::leftright(int w, int t, bool i, int x) : basic(w, t, i, x), onleft(true) {} 
	
bool leftright::bump_level() { 
    if (onleft) { 
        onleft=false;
        return true;
    } else {
        return false; 
    }
}

void leftright::set(int a, int t) { 
	row=a;
    if (onleft) {
        left=t-width;
        right=t-exclude;
    } else {
        left=t+exclude+1;
        right=t+width+1;
    } 
	restrain(); 
}	

/* The allaround class identifies all bin pairs in a surrounding square centred at the
 * bin pair of interest and excluding it. Here, bumping switches between the left and 
 * right sides if the stretch of bins is interrupted by the bin pairs to be excluded.
 */

allaround::allaround(int w, int t, bool i, int x) : basic(w, t, i, x), rangemode(FULLRANGE) { level=-w; } 

bool allaround::bump_level () { 
    if (rangemode==FULLRANGE) {
        if (level==-exclude-1) { 
            rangemode=LEFTSIDE;
        } else if (level >= width) { 
            return false; 
        }
        ++level;
    } else if (rangemode==LEFTSIDE) {
        rangemode=RIGHTSIDE; 
    } else if (rangemode==RIGHTSIDE) {
        if (level==exclude) {
            rangemode=FULLRANGE;
        } else {
            rangemode=LEFTSIDE;
        }
        ++level;
    }
    return true;
}

void allaround::set(int a, int t) {
	row=a+level;
    if (rangemode==FULLRANGE) {
        left=t-width;
        right=t+width+1;
    } else if (rangemode==LEFTSIDE) {
        left=t-width;
        right=t-exclude;
    } else if (rangemode==RIGHTSIDE) {
    	left=t+exclude+1;
        right=t+width+1;
    } 
    restrain();		
}

