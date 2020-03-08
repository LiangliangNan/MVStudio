
#ifndef _BASIC_ASSERTION_H_
#define _BASIC_ASSERTION_H_

#include "basic_common.h"
#include <string>



void BASIC_API ogf_assertion_failed(
									const std::string& condition_string,
									const std::string& file, int line ) ;

void BASIC_API ogf_range_assertion_failed(
	double value, double min_value, double max_value, 
	const std::string& file, int line ) ;

void BASIC_API ogf_should_not_have_reached(
	const std::string& file, int line ) ;


//________________________________________________________________________________


// Three levels of assert:
// use ogf_assert() and ogf_range_assert()               for non-expensive asserts
// use ogf_debug_assert() and ogf_debug_range_assert()   for expensive asserts
// use ogf_parano_assert() and ogf_parano_range_assert() for very exensive asserts

#define ogf_assert(x) {									\
	if(!(x)) {											\
	::ogf_assertion_failed(#x,__FILE__, __LINE__) ;		\
	}													\
} 

#define ogf_range_assert(x,min_val,max_val) {			\
	if(((x) < (min_val)) || ((x) > (max_val))) {		\
	::ogf_range_assertion_failed(x, min_val, max_val,	\
	__FILE__, __LINE__									\
	) ;													\
	}													\
}

#define ogf_assert_not_reached {						\
	::ogf_should_not_have_reached(__FILE__, __LINE__) ;	\
}

#ifndef NDEBUG
#define ogf_debug_assert(x) ogf_assert(x)
#define ogf_debug_range_assert(x,min_val,max_val) ogf_range_assert(x,min_val,max_val)
#else
#define ogf_debug_assert(x) 
#define ogf_debug_range_assert(x,min_val,max_val) 
#endif

#ifdef	PARANOID_DEBUG
#define ogf_parano_assert(x) ogf_assert(x)
#define ogf_parano_range_assert(x,min_val,max_val) ogf_range_assert(x,min_val,max_val)
#else
#define ogf_parano_assert(x) 
#define ogf_parano_range_assert(x,min_val,max_val) 
#endif

#endif
