#ifndef _BASIC_CANVAS_H_
#define _BASIC_CANVAS_H_

#include <string>


class PointSet;

class Canvas 
{
public:
	Canvas() {}
	virtual ~Canvas() {}

	virtual void fit() = 0; // fit screen
	virtual void update_graphics() = 0;
	virtual void update_all() = 0;

	// the active object
	virtual	PointSet*	pointSet() const = 0;
};


#endif