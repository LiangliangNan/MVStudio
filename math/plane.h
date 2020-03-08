
#ifndef _MATH_GEOMETRY_PLANE__
#define _MATH_GEOMETRY_PLANE__


#include "vecg.h"


// A 3D Plane of equation a.x + b.y + c.z + d = 0
template <class FT> 
class GenericPlane3 
{
public:
	typedef vecng<3, FT>		Point ;
	typedef vecng<3, FT>		Vector ;
	typedef GenericLine<3, FT>	Line ;
	typedef GenericPlane3<FT>	thisclass ;

public:
	GenericPlane3(const Point& p1, const Point& p2, const Point& p3) ;
	GenericPlane3(const Point& p, const Vector& n) ;
	GenericPlane3(double a, double b, double c, double d) : a_(a), b_(b), c_(c), d_(d) {}
	GenericPlane3() {}

	Vector normal() const ;

	// return a point lying on this plane
	Point  point() const;

	FT  value(const Point& p) const { return (a_*p.x + b_*p.y + c_*p.z + d_) ; }

	// return values:
	//   POSITIVE: p is on the positive side
	//   NEGATIVE: p is on the negative side
	//   ZERO:	   p is belonging to the plane
	Sign  orient(const Point& p) const ;

	// the projection of p on this plane
	Point projection(const Point &p) const ;

	FT	  squared_ditance(const Point &p) const ;

	bool  intersection(const Line& line, Point& p) const ;

private:
	FT a_;
	FT b_;
	FT c_;
	FT d_;
} ;


template <class FT> inline
GenericPlane3<FT>::GenericPlane3(const Point& p1, const Point& p2, const Point& p3) {
	Vector n = cross(p2-p1, p3-p1) ;
	a_ = n.x ;
	b_ = n.y ;
	c_ = n.z ;
	d_ = -( a_*p1.x + b_*p1.y + c_*p1.z ) ;

#ifndef NDEBUG // degenerate case
	if (length(n) < 1e-15) {
		std::cerr << "degenerate plane constructed from 3 points:" << std::endl
			<< "\t(" << p1 << ")" << std::endl 
			<< "\t(" << p2 << ")" << std::endl
			<< "\t(" << p3 << ")" << std::endl;
	}
#endif
}

template <class FT> inline
GenericPlane3<FT>::GenericPlane3(const Point& p, const Vector& n) {
	Vector nn = normalize(n);
	a_ = nn.x ;
	b_ = nn.y ;
	c_ = nn.z ;
	d_ = -( a_*p.x + b_*p.y + c_*p.z ) ;

#ifndef NDEBUG // degenerate case
	if (length(n) < 1e-15) {
		std::cerr << "degenerate plane constructed from point (" 
			<< p << ") and normal (" << n << ")" << std::endl;
	}
#endif
}


template <class FT> inline
typename GenericPlane3<FT>::Vector GenericPlane3<FT>::normal() const { 
	Vector n = normalize(Vector(a_, b_, c_)); 

#ifndef NDEBUG // degenerate case
	if (length(n) < 1e-15) {
		std::cerr << "degenerate plane with normal: (" << n << ")" << std::endl;
	}
#endif
	return n;
}


template < class FT> inline
typename GenericPlane3<FT>::Point GenericPlane3<FT>::point() const {
	Point p(0, 0, 0);
	if (a_ != 0)		p.x = -d_ / a_;
	else if (b_ != 0)	p.y = -d_ / b_;
	else				p.z = -d_ / c_;
	return p;
}


// return values:
//   1: p is on the positive side
//  -1: p is on the negative side
//   0: the point p is and 0 if the point is belonging the plane.
template <class FT> inline
Sign GenericPlane3<FT>::orient(const Point& p) const {
	FT v = value(p);
	if(ogf_abs(v) < 1e-15)
		return ZERO ;

	return (v > 0.0 ? POSITIVE : NEGATIVE) ;
}


template <class FT> inline
typename GenericPlane3<FT>::Point GenericPlane3<FT>::projection(const Point& p) const {
	// the equation of the plane is Ax+By+Cz+D=0
	// the normal direction is (A,B,C)
	// the projected point is p-lambda(A,B,C) where
	// A(x-lambda*A) + B(y-lambda*B) + C(z-lambda*C) + D = 0

	FT num = a_*p.x + b_*p.y + c_*p.z + d_;
	FT den = a_*a_ + b_*b_ + c_*c_;
	FT lambda = num / den;

	FT x = p.x - lambda * a_;
	FT y = p.y - lambda * b_;
	FT z = p.z - lambda * c_;
	return Point(x, y, z);
}


template <class FT> inline
FT GenericPlane3<FT>::squared_ditance(const Point& p) const {
	FT v = value(p);
	return (v * v) / (a_ * a_ + b_ * b_ + c_ * c_);
}


template <class FT> inline
bool GenericPlane3<FT>::intersection(const Line& line, Point& p) const {
	Vector dir = line.direction();
	FT c = dot(dir, normal());
	if (ogf_abs(c) < 1e-15)
		return false;

	Point p0 = line.point();
	// p = p0 + dir * t
	// equation: p is in the plane (so we first compute t)
	FT t = - (a_ * p0.x + b_ * p0.y + c_ * p0.z + d_) / (a_ * dir.x + b_ * dir.y + c_ * dir.z);
	p = p0 + dir * t;
	return true;
}

#endif
