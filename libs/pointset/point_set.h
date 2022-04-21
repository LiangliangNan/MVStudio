#ifndef _POINT_SET_H_
#define _POINT_SET_H_



#include "../basic/basic_types.h"
#include "../basic/object.h"
#include "../math/math_types.h"

#include <list>

class PointSet : public Object
{
public:
	PointSet();
	~PointSet();

	static std::string title() { return "PointSet"; }

	std::size_t  num_points() const { return points_.size(); }

	std::vector<vec3f>& points() { return points_; }
	std::vector<vec3f>& colors() { return colors_; }
	std::vector<vec3f>& normals() { return normals_; }
	const std::vector<vec3f>& points() const { return points_; }
	const std::vector<vec3f>& colors() const { return colors_; }
	const std::vector<vec3f>& normals() const { return normals_; }

	bool    has_normals() const { return normals_.size() > 0 && normals_.size() == points_.size(); }
	bool	has_colors() const  { return colors_.size() > 0 && colors_.size() == points_.size(); }

	void	delete_points(const std::vector<unsigned int>& indices);

	Box3f	bounding_box() const;

private:
	std::vector<vec3f>  points_;
	std::vector<vec3f>  colors_;
	std::vector<vec3f>  normals_;

};


#endif // _POINT_SET_H_
