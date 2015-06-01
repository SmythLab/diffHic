#include "neighbours.h"

basic::basic(int w, int t, bool i) : level(0), width(w), tlen(t), intra(i), remove_self(false) {}

bool basic::discard_self() { return (level==0 && remove_self); }

void basic::restrain () {
	if (left < 0) { left=0; }
	if (intra) {
		if (right > row) { 
			right=row+1; 
			if (left > right) { left=right; }
		}
	} else if (right > tlen) { right=tlen; } // For intra's, right will hit diagonal; no need to worry about tlen.
}

bottomright::bottomright(int w, int t, bool i) : basic(w, t, i) { level=-w; }

bool bottomright::bump_level () { 
	if (level >= 0) { return false; }
	++level;
	return true;
}

void bottomright::set(int a, int t) {
	row=a+level;
	left=(level ? t : t+1);
	right=t+width+1; 
	restrain();
}

updown::updown(int w, int t, bool i) : basic(w, t, i) { level=-w; }
	
bool updown::bump_level() { 
	if (level >= width) { return false; }
	++level;
	if (level==0) { ++level; }
	return true; 
}

void updown::set(int a, int t) {
	row=a+level;
	left=t;
	right=(level==0 ? t : t+1); 
	restrain();
}

leftright::leftright(int w, int t, bool i) : basic(w, t, i) { remove_self=true; }
	
bool leftright::bump_level() { return false; }

void leftright::set(int a, int t) { 
	row=a;
	left=t-width;
	right=t+width+1; 
	restrain(); 
}	

allaround::allaround(int w, int t, bool i) : basic(w, t, i) { 
	remove_self=true; 
	level=-w;
}

bool allaround::bump_level () { 
	if (level >= width) { return false; }
	++level;
	return true;
}

void allaround::set(int a, int t) {
	row=a+level;
	left=t-width;
	right=t+width+1;
	restrain();		
}

