
#ifndef _BASIC_TYPES_H_
#define _BASIC_TYPES_H_

#include "basic_common.h"
#include <string>
#include <vector>

#include <cmath>
#include <sstream>

#if (defined _WIN32) || (defined _WIN64)

#include <windows.h>
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#else

#include <unistd.h>

#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



//____________________________________________________________________________


/*
* A Namespace gathering typedefs for memory management.
*/

namespace Memory {

	typedef unsigned char  byte ;
	typedef unsigned char  word8 ;
	typedef unsigned short word16 ;
	typedef unsigned int   word32 ;

	typedef byte* pointer ;

	inline void clear(void* addr, size_t size) {
		::memset(addr, 0, size) ;
	}

	inline void copy(void* to, const void* from, size_t size) {
		::memcpy(to, from, size) ;
	}

} 

//_______________________________________________________________________

/**
* A namespace gathering typedefs
* corresponding to numbers. These types
* names have the form (u)int<size> or float<size>,
* where the (optional) u denotes an unsigned type,
* and the size is in bits.
*/

namespace Numeric {

	typedef char					int8 ;
	typedef short					int16 ;
	typedef int						int32 ;
#ifdef WIN32
	typedef __int64					int64 ;
#else
	typedef long long int			int64 ;	
#endif

	typedef unsigned char			uint8 ;
	typedef unsigned short			uint16 ;
	typedef unsigned int			uint32 ;

#ifdef WIN32
	typedef unsigned __int64		uint64 ;
#else
	typedef unsigned long long int	uint64 ;
#endif    

	typedef float					float32 ;
	typedef double					float64 ;

	typedef void*					pointer;

	// -----------------------------------------------------------------

	extern BASIC_API const float big_float ;
	extern BASIC_API const float small_float ;
	extern BASIC_API const double big_double ;
	extern BASIC_API const double small_double ;

	bool BASIC_API is_nan(float32 x) ;
	bool BASIC_API is_nan(float64 x) ;

	int32	BASIC_API random_int32() ;
	float32	BASIC_API random_float32() ;
	float64	BASIC_API random_float64() ;
} 

//_______________________________________________________________________

namespace String {

	void BASIC_API split_string(
		const std::string& in, 
		char separator,
		std::vector<std::string>& out,
		bool skip_empty_fields = true ) ;

	void BASIC_API join_strings(
		const std::vector<std::string>& in,
		char separator,
		std::string& out
		) ;

	void BASIC_API join_strings(
		const std::vector<std::string>& in,
		const std::string& separator,
		std::string& out
		) ;


	std::string BASIC_API join_strings(
		const std::vector<std::string>& in,
		char separator
		) ;

	std::string BASIC_API join_strings(
		const std::vector<std::string>& in,
		const std::string& separator
		) ;

	// return the number of substrings that have been replaced.
	// NOTE: if 'isolated' set to true, only isolated (by spaces) 
	//       substrings will be replaced.
	int BASIC_API replace_substring(
		std::string& in, 
		const std::string& original_substring, 
		const std::string& new_substring,
		bool isolated = true
		) ;

	// -----------------------------------------------------------------

	template <typename T>
	std::string to_string(T v)	{ 
		std::ostringstream ss;
		ss << v << '\0' ;
		//ss << v;	
		return ss.str(); 
	}

	template <typename T>
	T to_number(const std::string& str) { 
		std::istringstream ss(str);
		T result;	
		return (ss >> result) ? result : 0; 
	}

	// -----------------------------------------------------------------

	void BASIC_API to_lowercase(std::string& in) ;
	void BASIC_API to_uppercase(std::string& in) ;

	inline std::string BASIC_API char_to_string(char c) {
		char s[2] ;
		s[0] = c ;
		s[1] = '\0' ;
		return std::string(s) ;
	}

	std::string BASIC_API quote(const std::string& s, char quotes = '\"') ;

	// format example: "Fri Jan 09 11:39:32 2015"
	std::string BASIC_API from_current_time() ;
} 

//_______________________________________________________________________

#define nil 0

//_______________________________________________________________________

template <class T> 
inline T ogf_max(T x1, T x2) {
	return x1 > x2 ? x1 : x2;
}

template <class T> 
inline T ogf_min(T x1, T x2) {
	return x1 < x2 ? x1 : x2;
}

template <class T>
inline T ogf_max(T x1, T x2, T x3) {
	return ogf_max(ogf_max(x1, x2), x3);
}

template <class T>
inline T ogf_min(T x1, T x2, T x3) {
	return ogf_max(ogf_min(x1, x2), x3);
}

enum Sign { NEGATIVE=-1, ZERO=0, POSITIVE=1 } ;

template <class T> 
inline Sign ogf_sgn(T x) {
	return (x > 0) ? POSITIVE : (
		(x < 0) ? NEGATIVE : ZERO
		);
}

template <class T> 
inline T ogf_abs(T x)  {
	return (x > 0) ? x : -x;
}

template <class T> 
inline T ogf_sqr(T x)  {
	return x*x;
}

template <class T> 
inline void ogf_clamp(T& x, T min, T max) {
	if(x < min) {
		x = min ;
	} else if(x > max) {
		x = max ;
	}
}


template <class T> 
inline void ogf_swap(T& x, T& y) {
	T z = x ;
	x = y ;
	y = z ;
}

template <class T> 
inline T clip_precision(T v, int p) {
	float tmp = std::pow(10.0f, p);
	return (static_cast<int>(v * tmp) / tmp);
}


/* Return the closest integer to x, rounding up */
template <class T>
int iround(T x) {
	if (x < 0.0) {
		return (int)(x - 0.5);
	}
	else {
		return (int)(x + 0.5);
	}
}


#endif
