
#ifndef _COUNTED_H_
#define _COUNTED_H_

#include "basic_common.h"
#include "basic_types.h"
#include "smart_pointer.h"
#include "assertions.h"



//____________________________________________________________________________

/**
* This is the base class to be used for objects having
* "reference count" memory management. They can be 
* referred to by using SmartPointer<T>, calling ref()
* and unref() when necessary.
* @see SmartPointer
*/

class BASIC_API Counted {

public:
	Counted() ;
	virtual ~Counted() ;

	void ref() const ;
	void unref() const ;
	bool is_shared() const ;

	static void ref(const Counted* counted) ;
	static void unref(const Counted* counted) ;

protected:
private:
	int nb_refs_ ;
} ;

//____________________________________________________________________________

inline Counted::Counted() : nb_refs_(0) {
}

inline void Counted::ref() const {
	Counted* non_const_this = (Counted *)this ;
	non_const_this->nb_refs_++ ;
}

inline void Counted::unref() const {
	Counted* non_const_this = (Counted *)this ;    
	non_const_this->nb_refs_-- ;

	ogf_assert(nb_refs_ >= 0) ;

	if(nb_refs_ == 0) {
		delete this ;
	}
}

inline bool Counted::is_shared() const {
	return (nb_refs_ > 1) ;
}



#endif
