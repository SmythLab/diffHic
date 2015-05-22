#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

struct basic {
	basic(int, int, bool);
	virtual void set (int, int)=0;
	virtual ~basic() {};
	virtual bool bump_level ()=0; 
	bool discard_self();
	int row, left, right;
protected:
	int level, width, tlen;
	bool intra, remove_self;	
	void restrain ();
};


struct bottomright : public basic { 
	bottomright(int, int, bool);
	~bottomright() {};
	bool bump_level();
	void set(int, int);
};

struct updown : public basic {
	updown(int, int, bool);
	~updown() {};
	bool bump_level();
	void set(int, int);
};

struct leftright : public basic {
	leftright(int, int, bool);
	~leftright() {};
	bool bump_level();
	void set(int, int);
};

struct allaround : public basic {
	allaround(int, int, bool);
	~allaround() {};
	bool bump_level();
	void set(int, int);
};

#endif
