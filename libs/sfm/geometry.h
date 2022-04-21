
/* Geometric primitives */

#ifndef _SFM_GEOMETRY_H_
#define _SFM_GEOMETRY_H_

#include "../math/math_types.h"

#include <vector>


namespace sfm {


	typedef std::pair<int, int>		ImageKey;
	typedef std::vector<ImageKey>	ImageKeyVector;

	/* Data for tracks */
	class TrackData {
	public:
		TrackData() : extra(-1) {}
		TrackData(ImageKeyVector vs) : views(vs), extra(-1) { }

		/* Read/write routines */
		void read(FILE *f);
		void write(FILE *f);

		ImageKeyVector	views;
		int				extra;
	};


	/* Data for 3D points */
	class PointData {
	public:
		PointData() { fixed = false; }

		/* Write coordinates*/
		void write(FILE *f);

		vec3d	pos;	/* 3D position of the point */
		vec3d	norm;	/* Estimated normal for this point */
		vec3i	color;	/* Color of the point */
		double	conf;	/* Confidence in this point */

		ImageKeyVector	views;	/* View / keys corresponding to this point */
		bool	fixed;			/* Should this point be fixed during bundle adjustment? */

		float*	desc;			/* Descriptor for this point */
		int		num_vis;		/*number of images that see this point*/
		int		ref_image;		/* Reference image */
	};

}

#endif /* _SFM_GEOMETRY_H_ */
