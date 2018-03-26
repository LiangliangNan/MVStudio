

#ifndef _MATH_GEOMETRY_SEGMENT_H_
#define _MATH_GEOMETRY_SEGMENT_H_

#include "math_common.h"
#include "vecg.h"
#include "line.h"


template <int DIM, class FT> 
class GenericSegment {
public:
	typedef vecng<DIM, FT>			Point ;
	typedef vecng<DIM, FT>			Vector ;
	typedef GenericLine<DIM, FT>	Line ;
	typedef GenericSegment<DIM, FT> thisclass ;

public:
    GenericSegment(const Point& s, const Point& t) ;
	GenericSegment()  {}

	const Point& source() const { return s_ ; }
	const Point& target() const { return t_ ; }
	void  set_source(const Point& s) { s_ = s; }
	void  set_target(const Point& t) { t_ = t; }

	Line  supporting_line() const { return Line::from_two_points(s_, t_); }

	Vector to_vector() const { return t_ - s_ ; }

	// the projection of p on the supporting line 
	Point  projection(const Point &p) const { 
		Vector dir = normalize(t_ - s_);
		return (s_ + dir * dot(p - s_, dir)) ;
	}

	// test if the projection of p is inside the two end points
	bool   projected_inside(const Point& p) const {
		return (dot(s_ - p, t_ - p) < 0);
	}

	FT squared_ditance(const Point &p) const { 
		if (projected_inside(p))
			return distance2(projection(p), p);
		else {
			FT ds = distance2(s_, p);
			FT dt = distance2(t_, p);
			return ogf_min(ds, dt);
		}
	}

private:
	Point s_ ;
	Point t_ ;
} ;



template<int DIM, class FT> inline
GenericSegment<DIM, FT>::GenericSegment(const Point& s, const Point& t) : s_(s), t_(t) {
#ifndef NDEBUG // degenerate case
    if (distance2(s, t) < 1e-15) {
        std::cerr << "degenerate segment constructed from 2 points:" << std::endl
            << "\t(" << s << ")" << std::endl
            << "\t(" << t << ")" << std::endl;
    }

#endif
}


#endif
