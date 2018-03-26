
#ifndef __MATH_LINEAR_ALGEBRA_VECTOR__
#define __MATH_LINEAR_ALGEBRA_VECTOR__

#include "../basic/arrays.h"
#include "../basic/basic_types.h"
#include "../basic/assertions.h"

#include <iomanip>
#include <iostream>



class Vector : public Array1d<double> {
public:
	typedef Array1d<double> superclass ;
public:
	void zero() { Memory::clear(superclass::data(), superclass::size() * sizeof(double)); } 
	Vector() { }
	Vector(unsigned int n, unsigned int alignment = 1) : superclass(n, alignment) { zero(); }
	Vector(const Vector& rhs) {
		allocate(rhs.size(), rhs.alignment()) ;
		Memory::copy(superclass::data(), rhs.data(), superclass::size() * sizeof(double)) ;
	}
	void allocate(unsigned int n, unsigned int alignment=1) { superclass::allocate(n, alignment); zero() ; } 
	void operator += (const Vector& rhs) {
		ogf_debug_assert(rhs.size() == superclass::size()) ;
		for(unsigned int i=0; i<superclass::size(); i++) {
			(*this)(i) += rhs(i) ;
		}
	}
	void operator -= (const Vector& rhs) {
		ogf_debug_assert(rhs.size() == superclass::size()) ;
		for(unsigned int i=0; i<superclass::size(); i++) {
			(*this)(i) -= rhs(i) ;
		}
	}
	Vector& operator=(const Vector& rhs) {
		ogf_debug_assert(rhs.size() == superclass::size()) ;
		Memory::copy(superclass::data(), rhs.data(), superclass::size() * sizeof(double)) ;
		return *this ;
	}

	// used by operator << and may be used to print C-Style vectors
	static void print(std::ostream& out, unsigned int n, const double* x) {
		out << "[ " ;
		for(unsigned int i=0; i<n; i++) {
			out << std::setw(8) << x[i] << " " ;
		}
		out << "]" ;
	}
} ;

namespace Numeric {
	bool has_nan(const Vector& v) ;
} 

inline std::ostream& operator<<(std::ostream& out, const Vector& v) {
	v.print(out, v.size(), v.data()) ;
	return out ;
}

#endif
