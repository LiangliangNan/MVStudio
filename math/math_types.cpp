#include "math_types.h"
#include "../basic/basic_types.h"
#include "../basic/assertions.h"



namespace Geom {

	double triangle_beauty(
		const vec2d& p1, const vec2d& p2, const vec2d& p3
		) {
			vec2d p1p2 = p2 - p1 ;
			vec2d p2p3 = p3 - p2 ;
			vec2d p3p1 = p1 - p3 ;
			double l12 = ::sqrt(dot(p1p2, p1p2)) ;
			double l23 = ::sqrt(dot(p2p3, p2p3)) ;
			double l31 = ::sqrt(dot(p3p1, p3p1)) ;
			double d = 0.5 * (l12 + l23 + l31) ;
			return 8.0 * (d - l12) * (d - l23) * (d - l31) / (l12 * l23 * l31) ;
	}

	double triangle_beauty(
		const vec3d& p1, const vec3d& p2, const vec3d& p3
		) {
			vec3d p1p2 = p2 - p1 ;
			vec3d p2p3 = p3 - p2 ;
			vec3d p3p1 = p1 - p3 ;
			double l12 = ::sqrt(dot(p1p2, p1p2)) ;
			double l23 = ::sqrt(dot(p2p3, p2p3)) ;
			double l31 = ::sqrt(dot(p3p1, p3p1)) ;
			double d = 0.5 * (l12 + l23 + l31) ;
			return 
				8.0 * (d - l12) * (d - l23) * (d - l31) / (l12 * l23 * l31) ;
	}


	static inline void perp(const vec2d& p1, const vec2d& p2, vec2d& p, vec2d& v) {
		v = p2 - p1 ;
		v = vec2d(-v.y, v.x) ;
		p = barycenter(p1,p2) ;
	}

	vec2d triangle_circum_center(
		const vec2d& q1, const vec2d& q2, const vec2d& q3
		) {
			vec2d p1,p2 ;
			vec2d v1,v2 ;
			perp(q1,q2,p1,v1) ;
			perp(q1,q3,p2,v2) ;
			return segments_intersection_pv(p1,v1,p2,v2) ;
	}  

	vec3d perpendicular(const vec3d& V) {
		int min_index = 0 ;
		double c = ::fabs(V[0]) ;
		double cur = ::fabs(V[1]) ;
		if(cur < c) {
			min_index = 1 ;
			c = cur ;
		}
		cur = ::fabs(V[2]) ;
		if(cur < c) {
			min_index = 2 ;
		}
		vec3d result ;
		switch(min_index) {
			case 0:
				result = vec3d(0, -V.z, V.y) ;
				break ;
			case 1:
				result = vec3d(V.z, 0, -V.x) ;
				break ;
			case 2:
				result = vec3d(-V.y, V.x, 0) ;
				break ;
		}
		return result ;
	}

	inline double det3x3(
		double a00,  double a01,  double a02,
		double a10,  double a11,  double a12,
		double a20,  double a21,  double a22
		) {
			double m01 = a00*a11 - a10*a01;
			double m02 = a00*a21 - a20*a01;
			double m12 = a10*a21 - a20*a11;
			return m01*a22 - m02*a12 + m12*a02;
	}

	vec3d tetra_circum_center(
		const vec3d& p, const vec3d& q, 
		const vec3d& r, const vec3d& s
		) {
			vec3d qp = q-p ;
			double qp2 = length2(qp) ;
			vec3d rp = r-p ;
			double rp2 = length2(rp) ;
			vec3d sp = s-p ;
			double sp2 = length2(sp) ;

			double num_x = det3x3(
				qp.y,qp.z,qp2,
				rp.y,rp.z,rp2,
				sp.y,sp.z,sp2
				);

			double num_y = det3x3(
				qp.x,qp.z,qp2,
				rp.x,rp.z,rp2,
				sp.x,sp.z,sp2
				) ;

			double num_z = det3x3(
				qp.x,qp.y,qp2,
				rp.x,rp.y,rp2,
				sp.x,sp.y,sp2
				) ;

			double den = det3x3(
				qp.x,qp.y,qp.z,
				rp.x,rp.y,rp.z,
				sp.x,sp.y,sp.z
				) ;

			ogf_assert(::fabs(den) > 1e-30) ;

			den *= 2.0 ;

			return vec3d(
				p.x + num_x / den,
				p.y - num_y / den,
				p.z + num_z / den 
				) ;

	}


}
