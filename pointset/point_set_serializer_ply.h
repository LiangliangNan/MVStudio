
#ifndef _POINT_SERIALIZER_PLY_H_
#define _POINT_SERIALIZER_PLY_H_

#include <string>

class PointSet;
class PointSetSerializer_ply
{
public:
	static std::string title() { return "[PointSetSerializer_ply]: "; }

	static PointSet* load(const std::string& file_name) ;
	static bool		 save(const std::string& file_name, const PointSet* pset) ;
} ;

#endif

