#ifndef _BASIC_OBJECT_H_
#define _BASIC_OBJECT_H_

#include "basic_common.h"
#include <string>


class Canvas;

class BASIC_API Object 
{
public:
	Object() ;
	virtual ~Object() ;

	/**
	* If the mesh is a Object, updates the graphics display. This is 
	* used to animate the algorithms (for educational purposes)
	**/
	void fit();
	void update_graphics();
	void update_all();

	Canvas* canvas() const { return canvas_ ; }
	void set_canvas(Canvas* cvs) { canvas_ = cvs; } 

	const std::string& name() const { return name_ ; }
	void set_name(const std::string& n) { name_ = n ; }

protected:
	Canvas*		canvas_;
	std::string		name_; 
};


#endif