#ifndef _glQuickText_H_
#define _glQuickText_H_

#include "opengl_common.h"


/*
*   Written 2004 by <mgix@mgix.com>
*   This code is in the public domain
*   See http://www.mgix.com/snippets/?GLQuickText for details
*/

/* example usage:
Bbox3d box = Geom::map_bbox(mesh_);
std::string text = "My object \n vertices: %d, \n facets: %d";
glQuickText::printfAt(box.xmin(), box.ymin(), box.zmin(), 1.0f, text.c_str(), mesh_->size_of_vertices(), mesh_->size_of_facets()) ;
*/

class OPENGL_API GLQuickText
{
public:

	static void stringBox(
		double      *box,
		double      scale,
		const char  *format,
		...
		);

	static void printfAt(
		double      xPos,
		double      yPos,
		double      zPos,
		double      scale,
		const char  *format,
		...
		);

	static double getFontHeight(
		double scale = 1.0
		);
};

#endif // _glQuickText_H_

