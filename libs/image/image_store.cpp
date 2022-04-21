
#include "image_store.h"


ImageStore::ImageStore(int bytes_per_pixel, int dim_x) {
	bytes_per_pixel_ = bytes_per_pixel ;
	size_[0] = dim_x ;
	size_[1] = 1 ;
	size_[2] = 1 ;
	dimension_ = 1 ;
	base_mem_ = new Memory::byte[bytes()] ;
	Memory::clear(base_mem_, bytes()) ;
}

ImageStore::ImageStore(int bytes_per_pixel, int dim_x, int dim_y) {
	bytes_per_pixel_ = bytes_per_pixel ;
	size_[0] = dim_x ;
	size_[1] = dim_y ;
	size_[2] = 1 ;
	dimension_ = 2 ;
	base_mem_ = new Memory::byte[bytes()] ;
	Memory::clear(base_mem_, bytes()) ;
}

ImageStore::ImageStore(int bytes_per_pixel, int dim_x, int dim_y, int dim_z) {
	bytes_per_pixel_ = bytes_per_pixel ;
	size_[0] = dim_x ;
	size_[1] = dim_y ;
	size_[2] = dim_z ;
	dimension_ = 3 ;
	base_mem_ = new Memory::byte[bytes()] ;
	Memory::clear(base_mem_, bytes()) ;
}

ImageStore::~ImageStore() {
	delete[] base_mem_ ;
}



