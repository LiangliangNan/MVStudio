
#ifndef _MATH_TYPES_VECG_DYNAMIC_H_
#define _MATH_TYPES_VECG_DYNAMIC_H_

#include "math_common.h"
#include "../basic/basic_types.h"
#include "../basic/assertions.h"
#include <iostream>
#include <cfloat>


// vecng - dynamic



template<int DIM, class T> 
class vecng;


template <class T> 
class vecngd {
public:
	typedef vecngd<T> thisclass ;

	vecngd():dim_(0),data_(nil) { 			
	}

	explicit vecngd(unsigned int dim):dim_(dim),data_(new T[dim]) { 
		for(unsigned int i=0; i<dim; i++) { data_[i] = T(0); } 
	}

	~vecngd() {
		delete[] data_;
		data_ = nil ;
	}

	void resize(unsigned int dim) {
		if (dim_!=dim) {				
			delete[] data_;
			dim_=dim;
			data_ = new T[dim_];
		}
	}

	template<class T2> vecngd(const vecngd<T2>& rhs):dim_(rhs.dimension()),data_(new T[rhs.dimension()]) {
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] = T(rhs[i]) ;
		}
	}

	explicit vecngd(const vecngd<T>& rhs):dim_(rhs.dimension()),data_(new T[rhs.dimension()]) {			
		memcpy(data_, rhs.data(), dim_*sizeof(T));
	}

	template<class T2> vecngd(const T2* rhs, unsigned int dim):dim_(dim),data_(new T[dim]) {
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] = T(rhs[i]) ;
		}
	}
	vecngd(const T* rhs, unsigned int dim):dim_(dim),data_(new T[dim]) {
		memcpy(data_, rhs, dim_*sizeof(T));
	}
	template<class T2> vecngd(const T2* rhs, unsigned int dim, unsigned int stride):dim_(dim),data_(new T[dim]) {
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] = T(rhs[i*stride]) ;
		}
	}

	template<int DIM, class T2> vecngd(const vecng<DIM, T2>& rhs):dim_(DIM),data_(new T[DIM]) {
		for(unsigned int i=0; i<(unsigned int)dim_; i++) {
			data_[i] = T(rhs[i]) ;
		}
	}
	template<int DIM> vecngd(const vecng<DIM, T>& rhs):dim_(DIM),data_(new T[DIM]) {
		memcpy(data_, rhs.data(), dim_*sizeof(T));
	}

	thisclass& operator=(const thisclass& rhs) {
		if(&rhs != this) {
			resize(rhs.dimension());
			memcpy(data_, rhs.data(), dim_*sizeof(T));
		}
		return *this ;
	}

	template<int DIM, class T2> thisclass& operator=(const vecng<DIM,T2>& rhs) {
		resize(DIM);
		for(unsigned int i=0; i<DIM; i++) {
			data_[i] = T(rhs.data()[i]);
		}
		return *this ;
	}

	T* data() { return data_ ; }

	const T* data() const { return data_ ; }

	unsigned int getDimension() const {return dim_;}; // kept for the moment, but deprecated (not "coding-convetions"-compliant)
	unsigned int dimension() const {return dim_;};

	inline T& operator[](unsigned int idx) {
		ogf_debug_assert(idx < dim_) ;
		return data()[idx] ;
	}

	inline const T& operator[](unsigned int idx) const {
		ogf_debug_assert(idx < dim_) ;
		return data()[idx] ;
	}

	inline T length2() const { 
		T result = T(0) ;
		for(unsigned int i=0; i<dim_; i++) {
			result += data_[i]*data_[i] ;
		}
		return result ;
	}

	inline T distance2(const thisclass &rhs) const { 
		ogf_debug_assert(rhs.getDimension() == dim_) ;
		T result = T(0) ;
		for(unsigned int i=0; i<dim_; i++) {
			T val = rhs.data_[i]-data_[i];
			result += val*val ;
		}
		return result ;
	}

	inline T length() const {
		return sqrt(length2()) ;
	}

	// operators
	inline thisclass& operator+=(const thisclass& v) { 
		ogf_debug_assert(v.dimension() == dim_) ;
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] += v.data_[i] ;
		}
		return *this ; 
	}

	inline thisclass& operator-=(const thisclass& v) { 
		ogf_debug_assert(v.dimension() == dim_) ;
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] -= v.data_[i] ;
		}
		return *this ; 
	}

	inline thisclass& operator*=(const thisclass& v) {
		ogf_debug_assert(v.dimension() == dim_) ;
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] *= v.data_[i] ;
		}
		return *this ; 
	}

	inline thisclass& operator/=(const thisclass& v) { 
		ogf_debug_assert(v.dimension() == dim_) ;
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] /= v.data_[i] ;
		}
		return *this ; 
	}

	template <class T2> inline thisclass& operator*=(T2 s) { 
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] *= T(s) ;
		}
		return *this ; 
	}

	template <class T2> inline thisclass& operator/=(T2 s) { 
		for(unsigned int i=0; i<dim_; i++) {
			data_[i] /= T(s) ;
		}
		return *this ; 
	}

	inline thisclass operator+ (const thisclass& v) const {
		ogf_debug_assert(v.dimension() == dim_) ;
		thisclass result(*this);
		for(unsigned int i=0; i<dim_; i++) {
			result.data_[i] += v.data_[i] ;
		}
		return result ;
	}

	inline thisclass operator- (const thisclass& v) const {
		ogf_debug_assert(v.dimension() == dim_) ;
		thisclass result(*this);
		for(unsigned int i=0; i<dim_; i++) {
			result.data_[i] -= v.data_[i] ;
		}
		return result ;
	}

	template <class T2> inline thisclass operator* (T2 s) const {			
		thisclass result(*this);
		for(unsigned int i=0; i<dim_; i++) {
			result.data_[i] *= T(s) ;
		}
		return result ;
	}

	template <class T2> inline thisclass operator/ (T2 s) const {			
		thisclass result(*this);
		for(unsigned int i=0; i<dim_; i++) {
			result.data_[i] /= T(s) ;
		}
		return result ;
	}

	inline thisclass operator- () const {
		thisclass result(dim_);
		for(unsigned int i=0; i<dim_; i++) {
			result.data_[i] = -data_[i] ;
		}
		return result ;
	}

private:
	T* data_;
	unsigned int dim_;
} ;


template <class T> inline T dot(const vecngd<T>& v1, const vecngd<T>& v2) {  
	ogf_debug_assert(v1.dimension() == v2.dimension()) ;
	T result = 0 ;
	for(unsigned int i=0; i<v1.dimension(); i++) {
		result += v1[i]*v2[i] ;
	}
	return result ;
}

template <class T> inline vecngd<T> operator-(const vecngd<T>& v1) { 
	vecngd<T> result(v1.dimension()) ;
	for(unsigned int i=0; i<v1.dimension(); i++) {
		result[i] = -v1[i] ;
	}
	return result ;
}

template <class T2, class T> inline vecngd<T> operator*(T2 s, const vecngd<T>& v) { 
	vecngd<T> result(v.dimension()) ;
	for(unsigned int i=0; i<v.dimension(); i++) {
		result[i] = T(s)*v[i] ;
	}
	return result ;
}

template <class T> inline vecngd<T> operator+(const vecngd<T>& v1, const vecngd<T>& v2) { 
	ogf_debug_assert(v1.dimension() == v2.dimension()) ;
	vecngd<T> result(v1.dimension()) ;
	for(unsigned int i=0; i<v1.dimension(); i++) {
		result[i] = v1[i]+v2[i] ;
	}
	return result ;
}

template <class T> inline vecngd<T> operator-(const vecngd<T>& v1, const vecngd<T>& v2) { 
	ogf_debug_assert(v1.dimension() == v2.dimension()) ;
	vecngd<T> result(v1.dimension()) ;
	for(unsigned int i=0; i<v1.dimension(); i++) {
		result[i] = v1[i]-v2[i] ;
	}
	return result ;
}

// Compatibility with GLSL
template <class T> inline T length(const vecngd<T>& v) { return v.length() ; }
template <class T> inline T length2(const vecngd<T>& v) { return v.length2() ; }
template <class T> inline T distance2(const vecngd<T>& v1, const vecngd<T>& v2) { 
	return v2.distance2(v1) ; 
}
template <class T> inline T distance(const vecngd<T>& v1, const vecngd<T>& v2) { 
	return length(v2 - v1) ; 
}
template <class T> inline vecngd<T> normalize(const vecngd<T>& v) { 
	return (T(1) / length(v)) * v ; 
}
template <class T> inline vecngd<T> mix(const vecngd<T>& v1, const vecngd<T>& v2, T s) { 
	return (T(1) - s) * v1 + s * v2 ; 
}
template <class T> inline bool operator<(const vecngd<T> &v1, const vecngd<T>& v2)
{
	ogf_debug_assert(v1.dimension() == v2.dimension()) ;
	for(unsigned int i=0; i<v1.dimension(); i++) {
		if(v1.data()[i] < v2.data()[i]) { return true ;  }
		if(v1.data()[i] > v2.data()[i]) { return false ; }
	}
	return false ;
}

template <class T> inline std::ostream& operator<<(std::ostream& out, const vecngd<T>& v) {
	out << v.dimension()<<" ";
	for (unsigned int i=0; i<v.dimension()-1; i++)
		out << v[i]<< "  ";
	if (v.dimension()>0)
		out << v[v.dimension()-1];
	return out;
}

template <class T> inline std::istream& operator>>(std::istream& in, vecngd<T>& v) {
	unsigned int dim;
	in>>dim;
	v.resize(dim);
	for (unsigned int i=0; i<dim; i++)
		in >> v[i] ;
	return in;
}


#endif

