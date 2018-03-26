
#ifndef _MATH_GEOMETRY_BOX_H_
#define _MATH_GEOMETRY_BOX_H_

#include "math_common.h"
#include "vecg.h"
#include "../basic/assertions.h"


template <class T> 
class GenericBox2d {
public:
	GenericBox2d() : initialized_(false), x_min_(1e30), y_min_(1e30), x_max_(-1e30), y_max_(-1e30) { }

	bool initialized() const { return initialized_ ; }
	void clear() { initialized_ = false ; }

	T x_min() const { ogf_assert(initialized_); return x_min_ ; }
	T y_min() const { ogf_assert(initialized_); return y_min_ ; }
	T x_max() const { ogf_assert(initialized_); return x_max_ ; }
	T y_max() const { ogf_assert(initialized_); return y_max_ ; }

	T min(unsigned axis) const { return (axis==0)?x_min_:y_min_; }
	T max(unsigned axis) const { return (axis==0)?x_max_:y_max_; }

	T width()  const { return x_max() - x_min() ; }
	T height() const { return y_max() - y_min() ; }

	vecng<2,T> center() const {
		return vecng<2,T>((x_max() + x_min())/2, (y_max() + y_min())/2) ;
	}
	double radius() const {
		return 0.5 * ::sqrt( ogf_sqr(x_max()-x_min()) + ogf_sqr(y_max()-y_min()) ) ;
	}

	void add_point(const vecng<2,T>& p) {
		if(!initialized_) {
			x_min_ = p.x ;
			y_min_ = p.y ;
			x_max_ = p.x ;
			y_max_ = p.y ;
			initialized_ = true ;
		} else {
			x_min_ = ogf_min(x_min_, p.x) ;
			y_min_ = ogf_min(y_min_, p.y) ;
			x_max_ = ogf_max(x_max_, p.x) ;
			y_max_ = ogf_max(y_max_, p.y) ;
		}
	}

	void add_box(const GenericBox2d<T>& b) {
		if(b.initialized()) {
			add_point(vecng<2,T>(b.x_min(), b.y_min())) ;
			add_point(vecng<2,T>(b.x_max(), b.y_max())) ;
		}
	}

private:
	bool initialized_ ;
	T x_min_ ;
	T y_min_ ;
	T x_max_ ;
	T y_max_ ;
} ;

//_________________________________________________________________________

template <class T>
class GenericBox3d {
public:
	GenericBox3d() : initialized_(false), 
		x_min_(1e30), y_min_(1e30), z_min_(1e30), 
		x_max_(-1e30), y_max_(-1e30), z_max_(-1e30) { 
	}

	bool initialized() const { return initialized_ ; }
	void clear() { initialized_ = false ; }

	T x_min() const { ogf_assert(initialized_); return x_min_ ; }
	T y_min() const { ogf_assert(initialized_); return y_min_ ; }
	T z_min() const { ogf_assert(initialized_); return z_min_ ; }
	T x_max() const { ogf_assert(initialized_); return x_max_ ; }
	T y_max() const { ogf_assert(initialized_); return y_max_ ; }
	T z_max() const { ogf_assert(initialized_); return z_max_ ; }

	T min(unsigned axis) const { return (axis==0)?x_min_:((axis==1)?y_min_:z_min_); }
	T max(unsigned axis) const { return (axis==0)?x_max_:((axis==1)?y_max_:z_max_); }

	T width()  const { return x_max() - x_min() ; }
	T height() const { return y_max() - y_min() ; }
	T depth()  const { return z_max() - z_min() ; }

	vecng<3,T> center() const {
		return vecng<3,T>(
			0.5*(x_max() + x_min()), 
			0.5*(y_max() + y_min()),
			0.5*(z_max() + z_min())
			) ;
	}

	double radius() const {
		return 0.5 * ::sqrt( 
			ogf_sqr(x_max()-x_min()) + ogf_sqr(y_max()-y_min()) + ogf_sqr(z_max()-z_min()) 
			) ;
	}

	void add_point(const vecng<3,T>& p) {
		if(!initialized_) {
			x_min_ = p.x ;
			y_min_ = p.y ;
			z_min_ = p.z ;
			x_max_ = p.x ;
			y_max_ = p.y ;
			z_max_ = p.z ;
			initialized_ = true ;
		} else {
			x_min_ = ogf_min(x_min_, p.x) ;
			y_min_ = ogf_min(y_min_, p.y) ;
			z_min_ = ogf_min(z_min_, p.z) ;
			x_max_ = ogf_max(x_max_, p.x) ;
			y_max_ = ogf_max(y_max_, p.y) ;
			z_max_ = ogf_max(z_max_, p.z) ;
		}
	}

	void add_box(const GenericBox3d<T>& b) {
		if(b.initialized()) {
			add_point(vecng<3,T>(b.x_min(), b.y_min(), b.z_min())) ;
			add_point(vecng<3,T>(b.x_max(), b.y_max(), b.z_max())) ;
		}
	}

private:
	bool initialized_ ;
	T x_min_ ;
	T y_min_ ;
	T z_min_ ;
	T x_max_ ;
	T y_max_ ;
	T z_max_ ;
} ;



#endif
