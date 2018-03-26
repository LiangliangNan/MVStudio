
#ifndef _IMAGE_IMAGE_STORE_H_
#define _IMAGE_IMAGE_STORE_H_

#include "image_common.h"
#include "../basic/basic_types.h"
#include "../basic/counted.h"



/**
* Internal storage associated with a Image. A ImageStore
* is a multi-dimensional array of bytes. Note that several
* Images can share the same ImageStore.
*/

class IMAGE_API ImageStore : public Counted {
public:
	ImageStore(int bytes_per_pixel, int dim_x) ;
	ImageStore(int bytes_per_pixel, int dim_x, int dim_y) ;
	ImageStore(int bytes_per_pixel, int dim_x, int dim_y, int dim_z) ;
	virtual ~ImageStore() ;

	/**
	* Returns the dimension of this image. An image can be 1D, 2D or 3D. 
	*/
	int dimension() const ;

	/**
	* @return the size of the image along the specified axis, or
	*   0 if axis is greater than the dimension of this image.
	* @param axis can be one of 0,1,2 for X,Y,Z respectively.
	*/
	int size(int axis) const ;

	/** 
	* number of bytes allocated by this ImageStore
	*/
	int bytes() const ;

	/** equivalent to size(0) */
	int width() const ;

	/** equivalent to size(1). Returns 1 for a 1D image. */
	int height() const ;

	/** equivalent to size(2). Returns 1 for a 1D or a 2D image. */
	int depth() const ;

	int bytes_per_pixel() const ;

	Memory::pointer base_mem() const ;

private:
	int bytes_per_pixel_ ;
	int dimension_ ;
	int size_[3] ;
	Memory::byte* base_mem_ ;
} ;

typedef SmartPointer<ImageStore> ImageStore_var ;

//_________________________________________________________


inline int ImageStore::dimension() const {
	return dimension_ ;
}

inline int ImageStore::size(int axis) const {
	ogf_assert(axis >= 0 && axis < 3) ;
	return size_[axis] ;
}

inline int ImageStore::width() const {
	return size_[0] ;
}

inline int ImageStore::height() const {
	return size_[1] ;
}

inline int ImageStore::depth() const {
	return size_[2] ;
}

inline int ImageStore::bytes_per_pixel() const {
	return bytes_per_pixel_ ;
}

inline Memory::pointer ImageStore::base_mem() const {
	return base_mem_ ;
}

inline int ImageStore::bytes() const {
	return size_[0] * size_[1] * size_[2] * bytes_per_pixel_ ;
}

#endif

