#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "diffhic.h"

struct basic {
	basic(int, int, bool, int);
	virtual void set (int, int)=0;
	virtual ~basic() {};
	virtual bool bump_level ()=0; 
	int row, left, right;
protected:
	int level, width, tlen, exclude;
	bool intra;
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

struct leftright1 : public basic {
	leftright1(int, int, bool, int);
	~leftright1() {};
	bool bump_level();
	void set(int, int);
};

struct leftright2 : public basic {
	leftright2(int, int, bool, int);
	~leftright2() {};
	bool bump_level();
	void set(int, int);
};

struct allaround1 : public basic {
	allaround1(int, int, bool, int);
	~allaround1() {};
	bool bump_level();
	void set(int, int);
};

struct allaround2 : public basic {
	allaround2(int, int, bool, int);
	~allaround2() {};
	bool bump_level();
	void set(int, int);
};

#endif
