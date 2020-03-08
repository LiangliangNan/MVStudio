
#ifndef _IMAGE_IMAGE_H_
#define _IMAGE_IMAGE_H_



#include "image_store.h"
#include "color.h"


class Image {
public:
	enum ColorEncoding {
		GRAY, INDEXED, RGB, BGR, RGBA,
		INT16, INT32, FLOAT32, FLOAT64,
		RGB_FLOAT32, RGBA_FLOAT32,
		YUV
	} ;

	Image() ;
	Image(ColorEncoding color_rep, int size_x) ;
	Image(ColorEncoding color_rep, int size_x, int size_y) ;
	Image(ColorEncoding color_rep, int size_x, int size_y, int size_z) ;
	Image(const Image* rhs) ;

	virtual void initialize(ColorEncoding color_rep, int size_x) ;
	virtual void initialize(ColorEncoding color_rep, int size_x, int size_y) ;
	virtual void initialize(ColorEncoding color_rep, int size_x, int size_y, int size_z) ;

	virtual ~Image() ;

	virtual void acquire() ;

	int dimension() const { return dimension_ ; }
	int size(int axis) const { 
		assert(axis >= 0 && axis < 3) ;
		return size_[axis] ;
	}
	int width() const  { return size_[0] ; }
	int height() const { return size_[1] ; }
	int depth() const  { return size_[2] ; }

	int bytes_per_pixel() const { return bytes_per_pixel_ ; }
	int nb_pixels() const { return size_[0] * size_[1] * size_[2] ; }
	int bytes() const { return nb_pixels() * bytes_per_pixel_ ; }

	ColorEncoding color_encoding() const { return color_encoding_ ; }
	void set_color_encoding(ColorEncoding x) { color_encoding_ = x ; }

	easy3d::Memory::pointer base_mem() const {
		return base_mem_ ;
	}

    easy3d::Memory::byte* base_mem_byte_ptr() const {
		return byte_ptr(base_mem_) ;
	}

    easy3d::Numeric::int16* base_mem_int16_ptr() const {
		return int16_ptr(base_mem_) ;
	}

    easy3d::Numeric::int32* base_mem_int32_ptr() const {
		return int32_ptr(base_mem_) ;
	}

    easy3d::Numeric::float32* base_mem_float32_ptr() const {
		return float32_ptr(base_mem_) ;
	}

    easy3d::Numeric::float64* base_mem_float64_ptr() const {
		return float64_ptr(base_mem_) ;
	}

	static int bytes_per_pixel_from_color_encoding(
		ColorEncoding rep
		) ;

    easy3d::Memory::pointer pixel_base(int x) {
		if (x < 0) x = 0;
		if (x >= size_[0]) x = size_[0] - 1;
		return base_mem() + x * factor_[0] ;
	}

    easy3d::Memory::pointer pixel_base(int x, int y) {
		if (x < 0) x = 0;
		if (y < 0) y = 0;
		if (x >= size_[0]) x = size_[0] - 1;
		if (y >= size_[1]) y = size_[1] - 1;
		return base_mem() + x * factor_[0] + y * factor_[1] ;
	}

    easy3d::Memory::pointer pixel_base(int x, int y, int z) {
		if (x < 0) x = 0;
		if (y < 0) y = 0;
		if (z < 0) z = 0;
		if (x >= size_[0]) x = size_[0] - 1;
		if (y >= size_[1]) y = size_[1] - 1;
		if (z >= size_[2]) z = size_[2] - 1;
		return base_mem() + x * factor_[0] + y * factor_[1] + z * factor_[2];
	}

	inline easy3d::Memory::byte* byte_ptr(easy3d::Memory::pointer ptr) const {
		assert(
			color_encoding_ == GRAY ||
			color_encoding_ == RGB ||
			color_encoding_ == BGR ||
			color_encoding_ == RGBA ||
			color_encoding_ == YUV
			) ;
		return ptr ;
	}

	inline easy3d::Numeric::int16* int16_ptr(easy3d::Memory::pointer ptr) const {
		assert(color_encoding_ == INT16) ;
		return (easy3d::Numeric::int16*)(ptr) ;
	}

	inline easy3d::Numeric::int32* int32_ptr(easy3d::Memory::pointer ptr) const {
		assert(color_encoding_ == INT32) ;
		return (easy3d::Numeric::int32*)(ptr) ;
	}

	inline easy3d::Numeric::float32* float32_ptr(easy3d::Memory::pointer ptr) const {
		assert(
			color_encoding_ == FLOAT32 ||
			color_encoding_ == RGB_FLOAT32 ||
			color_encoding_ == RGBA_FLOAT32
			) ;
		return (easy3d::Numeric::float32*)(ptr) ;
	}

	inline easy3d::Numeric::float64* float64_ptr(easy3d::Memory::pointer ptr) const {
		assert(color_encoding_ == FLOAT64) ;
		return (easy3d::Numeric::float64*)(ptr) ;
	}

protected:
	ColorEncoding	color_encoding_ ;
	ImageStore*	store_ ;

	int factor_[3] ;
    easy3d::Memory::pointer base_mem_ ;
	int dimension_ ;
	int size_[3] ;
	int bytes_per_pixel_ ;

private:
	Image(const Image& rhs) ;
	Image& operator=(const Image& rhs) ;
} ;


// function provided to flip an image vertically
void   flip_image(Image* img);

// Interpolate the color for the point (x, y) in the given image
Colorf pixel_interpolate(Image *img, double x, double y);



#endif
