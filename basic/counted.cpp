
#include "counted.h"



//____________________________________________________________________________

Counted::~Counted() {
	ogf_assert(nb_refs_ == 0) ;
}

void Counted::ref(const Counted* counted) {
	if(counted != nil) {
		counted->ref() ;
	}
}

void Counted::unref(const Counted* counted) {
	if(counted != nil) {
		counted->unref() ;
	}
}

