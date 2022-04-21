#include "point_set.h"
#include "../basic/logger.h"


PointSet::PointSet()
{
}

PointSet::~PointSet() {
}


Box3f PointSet::bounding_box() const {
	Box3f result;
	for (int i = 0; i < points_.size(); ++i) {
		result.add_point(points_[i]);
	}
	return result;
}