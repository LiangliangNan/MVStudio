
#include "image.h"
#include <cmath>

using namespace easy3d;

Image::Image() {
	bytes_per_pixel_ = 0 ;
	dimension_ = 0 ;
	size_[0] = 0 ;
	size_[1] = 0 ;
	size_[2] = 0 ;
	factor_[0] = 0 ;
	factor_[1] = 0 ;
	factor_[2] = 0 ;
	base_mem_ = nullptr ;
}

Image::Image(ColorEncoding color_rep, int size_x) {
	initialize(color_rep, size_x) ;
}

Image::Image(ColorEncoding color_rep, int size_x, int size_y) {
	initialize(color_rep, size_x, size_y) ;
}

Image::Image(ColorEncoding color_rep, int size_x, int size_y, int size_z) {
	initialize(color_rep, size_x, size_y, size_z) ;
}

Image::Image(const Image* rhs) {
	initialize(
		rhs->color_encoding(), rhs->width(), rhs->height(), rhs->depth()
		) ;
	::memcpy(
		base_mem(), rhs->base_mem(),
		bytes_per_pixel() * width() * height() * depth() 
		) ;
}

void Image::initialize(ColorEncoding color_rep, int size_x) {
	color_encoding_ = color_rep ;
	bytes_per_pixel_ = bytes_per_pixel_from_color_encoding(color_rep) ;
	dimension_ = 1 ;
	size_[0] = size_x ;
	size_[1] = 1 ;
	size_[2] = 1 ;
	factor_[0] = bytes_per_pixel_ ;
	factor_[1] = 0 ;
	factor_[2] = 0 ;
	store_ = new ImageStore(bytes_per_pixel_, size_x) ;
	base_mem_ = store_->base_mem() ;
}

void Image::initialize(ColorEncoding color_rep, int size_x, int size_y) {
	color_encoding_ = color_rep ;
	bytes_per_pixel_ = bytes_per_pixel_from_color_encoding(color_rep) ;
	size_[0] = size_x ;
	size_[1] = size_y ;
	size_[2] = 1 ;
	dimension_ = 2 ;
	if(size_y == 1) {
		dimension_ = 1 ;
	}
	factor_[0] = bytes_per_pixel_ ;
	factor_[1] = factor_[0] * size_x ;
	factor_[2] = 0 ;
	store_ = new ImageStore(bytes_per_pixel_, size_x, size_y) ;
	base_mem_ = store_->base_mem() ;
}

void Image::initialize(ColorEncoding color_rep, int size_x, int size_y, int size_z) {
	color_encoding_ = color_rep ;
	bytes_per_pixel_ = bytes_per_pixel_from_color_encoding(color_rep) ;
	size_[0] = size_x ;
	size_[1] = size_y ;
	size_[2] = size_z ;
	dimension_ = 3 ;
	if(size_z == 1) {
		dimension_ = 2 ;
		if(size_y == 1) {
			dimension_ = 1 ;
		}
	}
	factor_[0] = bytes_per_pixel_ ;
	factor_[1] = factor_[0] * size_x ;
	factor_[2] = factor_[1] * size_y ;
	store_ = new ImageStore(bytes_per_pixel_, size_x, size_y, size_z) ;
	base_mem_ = store_->base_mem() ;
}

void Image::acquire() {
}

Image::~Image() {
    delete store_;
}

int Image::bytes_per_pixel_from_color_encoding(ColorEncoding rep) {
	int result = 0 ;
	switch(rep) {
		case GRAY:
			result = 1 ;
			break ;
		case INDEXED: 
			result = 1 ;
			break ;
		case RGB:
		case BGR:
		case YUV:
			result = 3 ;
			break ;
		case RGBA: 
			result = 4 ;
			break ;
		case INT16:
			result = 2 ;
			break ;
		case INT32:
			result = 4 ;
			break ;
		case FLOAT32:
			result = 4 ;
			break ;
		case FLOAT64:
			result = 8 ;
			break ;
		case RGB_FLOAT32:
			result = 12 ;
			break ;
		case RGBA_FLOAT32:
			result = 16 ;
			break ;
	}

	return result ;
}


void flip_image(Image* image) {
	int bpp = image->bytes_per_pixel();
	int h = image->height();
	int w = image->width();
	int row_len = w * bpp;

	for (int j = 0; j < h / 2; j++) {
		// get a pointer to the two lines we will swap
        easy3d::Memory::pointer row1 = image->base_mem() + j * row_len;
		easy3d::Memory::pointer row2 = image->base_mem() + (h - 1 - j) * row_len;
		// for each point on line, swap all the channels
		for (int i = 0; i < w; i++) {
			for (int k = 0; k < bpp; k++) {
				ogf_swap(row1[bpp*i + k], row2[bpp*i + k]);
			}
		}
	}
}



 #define LERP(x0, x1, f0, f1, f2, f3) ((1.0 - (x1)) * ((1.0 - (x0)) * (f0) + (x0) * (f1)) + (x1) * ((1.0 - (x0)) * (f2) + (x0) * (f3)))
 
 // Interpolate the color for the point (x, y) in the given image
 Colorf pixel_interpolate(Image *img, double x, double y) {
 	int xf = (int)floor(x), yf = (int)floor(y);
 	double xp = x - xf, yp = y - yf;
 
 	Memory::pointer pixels[4];
 	pixels[0] = img->pixel_base(xf, yf);
 	pixels[1] = img->pixel_base(xf + 1, yf);
 	pixels[2] = img->pixel_base(xf, yf + 1);
 	pixels[3] = img->pixel_base(xf + 1, yf + 1);
 
 	double rd = LERP(xp, yp, *(pixels[0]), *(pixels[1]), *(pixels[2]), *(pixels[3]));
 	double gd = LERP(xp, yp, *(pixels[0] + 1), *(pixels[1] + 1), *(pixels[2] + 1), *(pixels[3] + 1));
 	double bd = LERP(xp, yp, *(pixels[0] + 2), *(pixels[1] + 2), *(pixels[2] + 2), *(pixels[3] + 2));
 
 	return Colorf((float)rd, (float)gd, (float)bd);
 }