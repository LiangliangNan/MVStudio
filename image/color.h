
#ifndef _COLOR_H_
#define _COLOR_H_

#include <cassert>
#include <iostream>

#include <easy3d/core/types.h>

template <class T> 
class GenericColor
{
public:
	GenericColor(T r=0, T g=0, T b=0, T a=1.0f) ;
	template<class T2> explicit GenericColor(const T2* v, T a=1.0f) {
		components_[0] = v[0] ;
		components_[1] = v[1] ;
		components_[2] = v[2] ;
		components_[3] = a ;
	}

	T r() const ;
	T g() const ;
	T b() const ;
	T a() const ;

	void set_r(T r) ;
	void set_g(T g) ;
	void set_b(T b) ;
	void set_a(T a) ;

	void set(T r, T g, T b, T a=1.0f) ; 	
	
	T& operator[](int i) ;
	const T& operator[](int i) const ;

	// Low-level access
	const T* data() const ;
	T*       data() ;

	// color key 
	bool operator<(const GenericColor& rhs) const {
		if(components_[0] < rhs.components_[0]) { return true ; }
		if(components_[0] > rhs.components_[0]) { return false ; }
		if(components_[1] < rhs.components_[1]) { return true ; }
		if(components_[1] > rhs.components_[1]) { return false ; }
		if(components_[2] < rhs.components_[2]) { return true ; }
		if(components_[2] > rhs.components_[2]) { return false ; }
		return (components_[3] < rhs.components_[3]) ; 
	}

private:
	T components_[4] ;
} ;



template <class T> inline
GenericColor<T>::GenericColor(T r, T g, T b, T a) {
	components_[0] = r ;
	components_[1] = g ;
	components_[2] = b ;
	components_[3] = a ;
}

template <class T> inline
T GenericColor<T>::r() const {
	return components_[0] ;
}

template <class T> inline
T GenericColor<T>::g() const {
	return components_[1] ;
}

template <class T> inline
T GenericColor<T>::b() const {
	return components_[2] ;
}

template <class T> inline
T GenericColor<T>::a() const {
	return components_[3] ;
}

// Low-level access
template <class T> inline
const T* GenericColor<T>::data() const { 
	return components_; 
}

template <class T> inline
T* GenericColor<T>::data() {
	return components_; 
}

template <class T> inline
void GenericColor<T>::set(T r, T g, T b, T a) {
	components_[0] = r ;
	components_[1] = g ;
	components_[2] = b ;
	components_[3] = a ;
}

template <class T> inline
void GenericColor<T>::set_r(T r) {
	components_[0] = r ;
}

template <class T> inline
void GenericColor<T>::set_g(T g) {
	components_[1] = g ;
}

template <class T> inline
void GenericColor<T>::set_b(T b) {
	components_[2] = b ;
}

template <class T> inline
void GenericColor<T>::set_a(T a) {
	components_[3] = a ;
}

template <class T> inline
T& GenericColor<T>::operator[](int i) {
	assert(i >= 0 && i <= 3) ;
	return components_[i] ;
}

template <class T> inline
const T& GenericColor<T>::operator[](int i) const {
	assert(i >= 0 && i <= 3) ;
	return components_[i] ;
}


template <class T> inline
std::ostream& operator<<(std::ostream& output, const GenericColor<T>& color) {
	return output << 
		color[0] << " " << color[1] << " " << color[2] << " " << color[3] ;
}

template <class T> inline
std::istream& operator>>(std::istream& input, GenericColor<T>& color) {
	return input >> color[0] >> color[1] >> color[2] >> color[3] ;
}



//_______________________ Colors

typedef GenericColor<easy3d::Numeric::int8>    Color_int8 ;
typedef GenericColor<easy3d::Numeric::uint8>   Color_uint8 ;
typedef GenericColor<easy3d::Numeric::int16>   Color_int16 ;
typedef GenericColor<easy3d::Numeric::uint16>  Color_uint16 ;
typedef GenericColor<easy3d::Numeric::int32>   Color_int32 ;
typedef GenericColor<easy3d::Numeric::uint32>  Color_uint32 ;
typedef GenericColor<easy3d::Numeric::float32> Color_float32 ;
typedef GenericColor<easy3d::Numeric::float64> Color_float64 ;

typedef		Color_float32	Colorf ;


#endif
