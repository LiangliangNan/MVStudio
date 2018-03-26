#include "vectord.h"



namespace Numeric {

	bool has_nan(const Vector& v) {
		for(unsigned int i=0; i<v.size(); i++) {
			if(is_nan(v(i))) {
				return true ;
			}
		}
		return false ;
	}

}


