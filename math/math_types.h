
#ifndef _MATH_TYPES_H_
#define _MATH_TYPES_H_


#include "vecg.h"
#include "box.h"
#include "line.h"
#include "segment.h"
#include "plane.h"
#include "matrix.h"

/**
 * Gathers different base types for geometric operations. 
 * Types defined here are points, vectors, lines, segments, 
 * boxes and matrices in 2D and 3D. 
 */

class GeometricTypes
{
public:

	//---------------------- Points in 2d

	typedef vecng<2, Numeric::int8>				Point2d_int8 ;
	typedef vecng<2, Numeric::int16>			Point2d_int16 ;
	typedef vecng<2, Numeric::int32>			Point2d_int32 ;
	typedef vecng<2, Numeric::float32>			Point2d_float32 ;
	typedef vecng<2, Numeric::float64>			Point2d_float64 ;

	//_______________________ Points in 3d

	typedef vecng<3, Numeric::int8>				Point3d_int8 ;
	typedef vecng<3, Numeric::int16>			Point3d_int16 ;
	typedef vecng<3, Numeric::int32>			Point3d_int32 ;
	typedef vecng<3, Numeric::float32>			Point3d_float32 ;
	typedef vecng<3, Numeric::float64>			Point3d_float64 ;
	
	//_______________________ Points in 4d

	typedef vecng<4, Numeric::int8>				Point4d_int8 ;
	typedef vecng<4, Numeric::int16>			Point4d_int16 ;
	typedef vecng<4, Numeric::int32>			Point4d_int32 ;
	typedef vecng<4, Numeric::float32>			Point4d_float32 ;
	typedef vecng<4, Numeric::float64>			Point4d_float64 ;

	//_______________________ Vectors in 2d

	typedef vecng<2, Numeric::int8>				Vector2d_int8 ;
	typedef vecng<2, Numeric::int16>			Vector2d_int16 ;
	typedef vecng<2, Numeric::int32>			Vector2d_int32 ;
	typedef vecng<2, Numeric::float32>			Vector2d_float32 ;
	typedef vecng<2, Numeric::float64>			Vector2d_float64 ;

	//_______________________ Vectors in 3d

	typedef vecng<3, Numeric::int8>				Vector3d_int8 ;
	typedef vecng<3, Numeric::int16>			Vector3d_int16 ;
	typedef vecng<3, Numeric::int32>			Vector3d_int32 ;
	typedef vecng<3, Numeric::float32>			Vector3d_float32 ;
	typedef vecng<3, Numeric::float64>			Vector3d_float64 ;
	
	//_______________________ Vectors in 4d

	typedef vecng<4, Numeric::int8>				Vector4d_int8 ;
	typedef vecng<4, Numeric::int16>			Vector4d_int16 ;
	typedef vecng<4, Numeric::int32>			Vector4d_int32 ;
	typedef vecng<4, Numeric::float32>			Vector4d_float32 ;
	typedef vecng<4, Numeric::float64>			Vector4d_float64 ;
	
	//_______________________ Lines in 2d

	typedef GenericLine<2, Numeric::int8>		Line2d_int8 ;
	typedef GenericLine<2, Numeric::int16>		Line2d_int16 ;
	typedef GenericLine<2, Numeric::int32>		Line2d_int32 ;
	typedef GenericLine<2, Numeric::float32>	Line2d_float32 ;
	typedef GenericLine<2, Numeric::float64>	Line2d_float64 ;
	
	//_______________________ Lines in 3d

	typedef GenericLine<3, Numeric::int8>		Line3d_int8 ;
	typedef GenericLine<3, Numeric::int16>		Line3d_int16 ;
	typedef GenericLine<3, Numeric::int32>		Line3d_int32 ;
	typedef GenericLine<3, Numeric::float32>	Line3d_float32 ;
	typedef GenericLine<3, Numeric::float64>	Line3d_float64 ;

	//_______________________ Segments in 2d

	typedef GenericSegment<2, Numeric::int8>    Segment2d_int8 ;
	typedef GenericSegment<2, Numeric::int16>   Segment2d_int16 ;
	typedef GenericSegment<2, Numeric::int32>   Segment2d_int32 ;
	typedef GenericSegment<2, Numeric::float32> Segment2d_float32 ;
	typedef GenericSegment<2, Numeric::float64> Segment2d_float64 ;

	//_______________________ Segments in 3d

	typedef GenericSegment<3, Numeric::int8>    Segment3d_int8 ;
	typedef GenericSegment<3, Numeric::int16>   Segment3d_int16 ;
	typedef GenericSegment<3, Numeric::int32>   Segment3d_int32 ;
	typedef GenericSegment<3, Numeric::float32> Segment3d_float32 ;
	typedef GenericSegment<3, Numeric::float64> Segment3d_float64 ;

	//__________________________ Boxes in 2d

	typedef GenericBox2d<Numeric::float32>		Box2d_float32 ;
	typedef GenericBox2d<Numeric::float64>		Box2d_float64 ;

	//__________________________ Boxes in 3d

	typedef GenericBox3d<Numeric::float32>		Box3d_float32 ;
	typedef GenericBox3d<Numeric::float64>		Box3d_float64 ;

	//_________________________ Planes in 3d

	typedef GenericPlane3<Numeric::float32>		Plane3d_float32 ;
	typedef GenericPlane3<Numeric::float64>		Plane3d_float64 ;

	//_________________________ Matrices 

	typedef Matrix<Numeric::float32, 2>			Matrix2_float32 ;
	typedef Matrix<Numeric::float64, 2>			Matrix2_float64 ;

	typedef Matrix<Numeric::float32, 3>			Matrix3_float32 ;
	typedef Matrix<Numeric::float64, 3>			Matrix3_float64 ;

	typedef Matrix<Numeric::float32, 4>			Matrix4_float32 ;
	typedef Matrix<Numeric::float64, 4>			Matrix4_float64 ;


} ;



//____________________ default types___________________

//---------------------- Point / vector

typedef vecng<2, Numeric::int32>			vec2i;
typedef vecng<3, Numeric::int32>			vec3i;

typedef vecng<2, Numeric::float32>			vec2f;
typedef vecng<3, Numeric::float32>			vec3f;

typedef vecng<2, Numeric::float64>			vec2d;
typedef vecng<3, Numeric::float64>			vec3d;

typedef Matrix<Numeric::float32, 2>			mat2f;
typedef Matrix<Numeric::float32, 3>			mat3f;
typedef Matrix<Numeric::float32, 4>			mat4f;

typedef Matrix<Numeric::float64, 2>			mat2d;
typedef Matrix<Numeric::float64, 3>			mat3d;
typedef Matrix<Numeric::float64, 4>			mat4d;

typedef GenericBox2d<Numeric::float32>		Box2f;
typedef GenericBox3d<Numeric::float32>		Box3f;

typedef GenericBox2d<Numeric::float64>		Box2d;
typedef GenericBox3d<Numeric::float64>		Box3d;


//_________________________________________________________


// This namespace gathers some global geometric utilities.
namespace Geom {


	/* Compute the mean of a set of vectors */
	inline vec3d vec_mean(int n, const vec3d *v) {
		vec3d mean(0.0, 0.0, 0.0);
		for (int i = 0; i < n; i++) {
			mean = mean + v[i];
		}
		return mean / n;
	}

	inline vec2d vec_mean(int n, const vec2d *v) {
		vec2d mean(0.0, 0.0);
		for (int i = 0; i < n; i++) {
			mean = mean + v[i];
		}
		return mean / n;
	}

    inline vec3d barycenter(const vec3d& p1, const vec3d& p2) {
        return vec3d(
            0.5f * (p1.x + p2.x),
            0.5f * (p1.y + p2.y),
            0.5f * (p1.z + p2.z)
        ) ;
    }

    inline vec2d barycenter(const vec2d& p1, const vec2d& p2) {
        return vec2d(
            0.5f * (p1.x + p2.x),
            0.5f * (p1.y + p2.y)
        ) ;
    }

    inline vec3d barycenter(
        const vec3d& p1, const vec3d& p2, const vec3d& p3
    ) {
        return vec3d(
            (p1.x + p2.x + p3.x) / 3.0f ,
            (p1.y + p2.y + p3.y) / 3.0f,
            (p1.z + p2.z + p3.z) / 3.0f
        ) ;
    }

    inline vec2d barycenter(
        const vec2d& p1, const vec2d& p2, const vec2d& p3
    ) {
        return vec2d(
            (p1.x + p2.x + p3.x) / 3.0f ,
            (p1.y + p2.y + p3.y) / 3.0f 
        ) ;
    }

    inline vec3d barycenter(
        const vec3d& p1, const vec3d& p2, const vec3d& p3, const vec3d& p4
    ) {
        return vec3d(
            0.25f * (p1.x + p2.x + p3.x + p4.x) ,
            0.25f * (p1.y + p2.y + p3.y + p4.y) ,
            0.25f * (p1.z + p2.z + p3.z + p4.z)
        ) ;
    }

#define DEG2RAD(d) ((d) * (M_PI / 180.0))
#define RAD2DEG(r) ((r) * (180.0 / M_PI))

    inline double cos_angle(const vec3d& a, const vec3d& b) {
        double na2 = length2(a) ;
        double nb2 = length2(b) ;
        return dot(a, b) / ::sqrt(na2 * nb2) ;
    }

    inline double angle(const vec3d& a, const vec3d& b) {
        return ::acos(cos_angle(a,b)) ;
    }


    inline double cos_angle(const vec2d& a, const vec2d& b) {
        double na2 = length2(a) ;
        double nb2 = length2(b) ;
        return dot(a, b) / ::sqrt(na2 * nb2) ;
    }

    inline double det(const vec2d& a, const vec2d& b) {
        return a.x*b.y - a.y*b.x ;            
    }

    
    /* returns the angle in the interval [-pi .. pi] */
    inline double angle(const vec2d& a, const vec2d& b) {
        return det(a,b) > 0         ? 
             ::acos(cos_angle(a,b)) : 
            -::acos(cos_angle(a,b)) ;
    }

    inline double triangle_area(
        const vec3d& p1, const vec3d& p2, const vec3d& p3
    ) {
        return 0.5 * length(cross(p2 - p1, p3 - p1)) ;
    }

    inline double triangle_signed_area(
        const vec2d& p1, const vec2d& p2, const vec2d& p3
    ) {
        return 0.5 * det(p2-p1,p3-p1) ;
    }

    inline double triangle_area(
        const vec2d& p1, const vec2d& p2, const vec2d& p3
    ) {
        return ::fabs(triangle_signed_area(p1,p2,p3)) ;
    }

    double triangle_beauty(
        const vec2d& p1, const vec2d& p2, const vec2d& p3
    ) ;

    double triangle_beauty(
        const vec3d& p1, const vec3d& p2, const vec3d& p3
    ) ;

    vec2d triangle_circum_center(
        const vec2d& p1, const vec2d& p2, const vec2d& p3
    ) ;

    inline bool exact_equal(const vec3d& a,const vec3d& b) {
        return a.x==b.x && a.y==b.y && a.z==b.z ;
    }

    inline bool exact_equal(const vec2d& a,const vec2d& b) {
        return a.x==b.x && a.y==b.y ;
    }

    inline bool is_nan(const vec3d& v) {
        return 
            Numeric::is_nan(v.x) ||
            Numeric::is_nan(v.y) ||
            Numeric::is_nan(v.z) ;
    }

    inline bool is_nan(const vec2d& v) {
        return 
            Numeric::is_nan(v.x) ||
            Numeric::is_nan(v.y) ;
    }

    inline double segments_intersection_parameter(
        const vec2d& p1, const vec2d& v1,
        const vec2d& p2, const vec2d& v2
    ) {
        double delta = v1.x * v2.y - v1.y * v2.x ;
        double t1 = (v2.y * (p2.x - p1.x) - v2.x * (p2.y - p1.y)) / delta ;
        return t1 ;
    }

    inline vec2d segments_intersection_pv(
        const vec2d& p1, const vec2d& v1,
        const vec2d& p2, const vec2d& v2
    ) {
        double t1 = segments_intersection_parameter(p1, v1, p2, v2) ;
        return p1 + t1 * v1 ;
    }

    inline vec2d segments_intersection_pp(
        const vec2d& p1, const vec2d& p2,
        const vec2d& p3, const vec2d& p4
    ) {
        return segments_intersection_pv(
            p1, (p2 -p1), p3, (p4-p3)
        ) ;
    }

    vec3d perpendicular(const vec3d& V) ;

    inline double tetra_signed_volume(
        const vec3d& p1, const vec3d& p2, 
        const vec3d& p3, const vec3d& p4
    ) {
        return dot(p2 - p1, cross(p3 - p1, p4 - p1)) / 6.0 ;
    }

    inline double tetra_volume(
        const vec3d& p1, const vec3d& p2, 
        const vec3d& p3, const vec3d& p4
    ) {
        return ::fabs(tetra_signed_volume(p1,p2,p3,p4)) ;
    }

    vec3d tetra_circum_center(
        const vec3d& p1, const vec3d& p2, 
        const vec3d& p3, const vec3d& p4
    ) ;

}


#endif


