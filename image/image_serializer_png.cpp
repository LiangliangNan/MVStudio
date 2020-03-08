
#include "image_serializer_png.h"
#include "image.h"
#include "png.h"
#include <easy3d/util/logging.h>

Image* ImageSerializer_png::serialize_read(const std::string& file_name) {
	FILE* in = fopen(file_name.c_str(), "rb") ;
	//FILE* in; fopen_s(&in, file_name.c_str(), "rb") ; // fopen_s doesn't exist under Linux and Mac
	if(in == nullptr) {
		LOG(ERROR) << "Could not open file: \'"
			<< file_name
			<< "\'" << std::endl ;
		return nullptr ;
	}

	png_structp png_ptr = png_create_read_struct(
		PNG_LIBPNG_VER_STRING,
		(png_voidp) NULL, (png_error_ptr) NULL, (png_error_ptr) NULL 
		) ;
	if ( png_ptr == nullptr ) {
		fclose(in) ;
		return nullptr ;
	} 

	png_infop info_ptr = png_create_info_struct ( png_ptr );
	if ( info_ptr == nullptr ) {
		png_destroy_read_struct (
			&png_ptr, (png_infopp)NULL, (png_infopp)NULL 
			) ;
		fclose(in) ;
		return nullptr ;
	}

	png_init_io ( png_ptr, in );

	// read header 
	int bit_depth, color_type, interlace_type ;
	png_uint_32 width, height;
	png_read_info ( png_ptr, info_ptr );
	png_get_IHDR ( 
		png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
		&interlace_type, NULL, NULL 
		) ;


	Image* result = nullptr;

	if(color_type == PNG_COLOR_TYPE_GRAY) {
		result = new Image(Image::GRAY, width, height) ;
	} else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
		result = new Image(Image::RGBA, width, height) ;
	} else {
		result = new Image(Image::RGB, width, height) ;
	}

	// If colormapped, convert to RGB
	// TODO: create a ColorMapped Image
	if ( 
		color_type == PNG_COLOR_TYPE_PALETTE || // expand colormapped to RGB
		color_type == PNG_COLOR_TYPE_GRAY       // we also expand less than 8 bpp grayscale images to 8 bpp.
		) {
			png_set_expand(png_ptr) ;
	}

	// Ignore alpha (for the moment)
	// TODO: read alpha if present
	//png_set_strip_alpha(png_ptr) ;

	// Read the image one line at a time
	png_bytep row_pointer = (png_bytep) png_malloc ( 
		png_ptr, png_get_rowbytes ( png_ptr, info_ptr ) 
		) ;
	row_pointer = (png_bytep) malloc ( result->width() * 4 );
	for (int row = height-1; row >= 0; row-- ) {
		png_read_rows ( png_ptr, &row_pointer, NULL, 1 );
		if(
			color_type == PNG_COLOR_TYPE_GRAY ||
			color_type == PNG_COLOR_TYPE_PALETTE
			) {
				memcpy ( 
					result->base_mem() + row * result->width(),
					row_pointer,  result->width() 
					) ;
		} else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
			memcpy ( 
				result->base_mem() + row * result->width() * 4,
				row_pointer,  result->width() * 4
				) ;
		} else {
			memcpy ( 
				result->base_mem() + row * result->width() * 3,
				row_pointer,  result->width() * 3  
				) ;
		}
	}
	free ( row_pointer );

	png_read_end ( png_ptr, info_ptr );
	png_destroy_read_struct ( &png_ptr, &info_ptr, (png_infopp)NULL );

	fclose(in) ;
	return result ;
}

bool ImageSerializer_png::read_supported() const {
	return true ;
}

bool ImageSerializer_png::serialize_write(
	const std::string& file_name, const Image* image
	) {
		if(
			image->color_encoding() != Image::RGB && 
			image->color_encoding() != Image::RGBA  &&
			image->color_encoding() != Image::GRAY
			) {
				LOG(ERROR)
					<< "PNG writing only supported for GRAY, RGB and RGBA color encoding"
					<< ", sorry" << std::endl ;
				return false ;
		}

		FILE* out = fopen(file_name.c_str(), "wb") ;
		//FILE* out; fopen_s(&out, file_name.c_str(), "wb") ; // fopen_s doesn't exist under Linux and Mac
		if(out == nullptr) {
			LOG(ERROR) << "Could not open file: \'"
				<< file_name
				<< "\'" << std::endl ;
			return false ;
		}

		png_structp png_ptr = png_create_write_struct(
			PNG_LIBPNG_VER_STRING,
			(png_voidp) NULL, (png_error_ptr) NULL, (png_error_ptr) NULL 
			) ;
		if ( png_ptr == nullptr ) {
			fclose(out) ;
			return false ;
		} 

		png_infop info_ptr = png_create_info_struct ( png_ptr );
		png_init_io( png_ptr, out );


		png_byte png_color_encoding ;

		switch(image->color_encoding()) {
		case Image::GRAY:
			png_color_encoding = PNG_COLOR_TYPE_GRAY ;
			break ;
		case Image::RGB : 
			png_color_encoding = PNG_COLOR_TYPE_RGB ;
			break ;
		case Image::RGBA : 
			png_color_encoding = PNG_COLOR_TYPE_RGBA ;
			break ;
		default:
			LOG(ERROR) << "should not reach here" ;
			break ;
		}

		png_set_IHDR ( 
			png_ptr, info_ptr, image->width(), image->height(),
			8, png_color_encoding, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE 
			);

		png_text comment;
		comment.compression = PNG_TEXT_COMPRESSION_NONE;
		comment.key = "Comment" ;
		comment.text = "Made with MeshStudio" ;
		png_set_text ( png_ptr, info_ptr, &comment, 1 );

		png_write_info ( png_ptr, info_ptr );
		for (int row = image->height()-1; row >= 0 ; row-- ) {
			png_write_row ( 
				png_ptr, 
				image->base_mem() + image->width() * row * image->bytes_per_pixel()
				);
		}
		png_write_end ( png_ptr, info_ptr );
		fflush ( out );
		png_destroy_write_struct ( &png_ptr, (png_infopp)NULL );
		fclose(out) ;

		return true ;
}

bool ImageSerializer_png::write_supported() const {
	return true ;
}

bool ImageSerializer_png::streams_supported() const {
	return false ;
}

bool ImageSerializer_png::binary() const {
	return true ;
}


bool ImageSerializer_png::query_image_size(const std::string& file_name, int& width, int& height) {
	FILE* in = fopen(file_name.c_str(), "rb");
	//FILE* in; fopen_s(&in, file_name.c_str(), "rb") ; // fopen_s doesn't exist under Linux and Mac
	if (in == nullptr) {
		LOG(ERROR) << "Could not open file: \'"
			<< file_name
			<< "\'" << std::endl;
		return false;
	}

	png_structp png_ptr = png_create_read_struct(
		PNG_LIBPNG_VER_STRING,
		(png_voidp)NULL, (png_error_ptr)NULL, (png_error_ptr)NULL
		);
	if (png_ptr == nullptr) {
		fclose(in);
		return false;
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == nullptr) {
		png_destroy_read_struct(
			&png_ptr, (png_infopp)NULL, (png_infopp)NULL
			);
		fclose(in);
		return false;
	}

	png_init_io(png_ptr, in);

	// read header 
	int bit_depth, color_type, interlace_type;
	png_uint_32 w, h;
	png_read_info(png_ptr, info_ptr);
	if (png_get_IHDR(
		png_ptr, info_ptr, &w, &h, &bit_depth, &color_type,
		&interlace_type, NULL, NULL
		))
		return true;
	else
		return false;
}