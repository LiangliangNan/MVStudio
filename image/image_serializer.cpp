
#include "image_serializer.h"
#include "image.h"
#include "../basic/logger.h"

#include <fstream>


Image* ImageSerializer::serialize_read(const std::string& file_name) {

	std::fstream::openmode mode = binary() ?
		(std::fstream::in | std::fstream::binary) :
	std::fstream::in ;

	std::ifstream input(file_name.c_str(), mode) ;
	if(!input) {
		Logger::err("ImageSerializer") 
			<< "could not open file\'" 
			<< file_name << "\'" << std::endl ;
		return nil ;
	}
	return serialize_read(input) ;
}

bool ImageSerializer::serialize_write(
									  const std::string& file_name, const Image* image
									  ) 
{
	std::fstream::openmode mode = binary() ?
	  (std::fstream::out | std::fstream::trunc | std::fstream::binary) :
	(std::fstream::out | std::fstream::trunc) ;

	std::ofstream output(file_name.c_str(), mode) ;

	if(!output) {
	  Logger::err("ImageSerializer") 
		  << "could not open file\'" 
		  << file_name << "\'" << std::endl ;
	  return false ;
	}

	return serialize_write(output, image) ;
}

Image* ImageSerializer::serialize_read(std::istream& stream) {
	bool implemented = false ;
	ogf_assert(implemented) ;
	return nil ;
}

bool ImageSerializer::serialize_write(
									  std::ostream& stream, const Image* image
									  )
{
	bool implemented = false ;
	ogf_assert(implemented) ;
	return false ;
}


bool ImageSerializer::binary() const {
	return true ;
}

bool ImageSerializer::streams_supported() const {
	return true ;
}

bool ImageSerializer::read_supported() const {
	return false ;
}

bool ImageSerializer::write_supported() const {
	return false ;
}

bool ImageSerializer::flip_image(Image& image) const {
	int bpp=image.bytes_per_pixel();
	int h = image.height() ;
	int w = image.width() ;
	int row_len = w * bpp ;

	for(int j=0; j< h/2; j++) {
		// get a pointer to the two lines we will swap
		Memory::pointer row1 = image.base_mem() + j * row_len ;
		Memory::pointer row2 = image.base_mem() + (h - 1 - j) * row_len ;
		// for each point on line, swap all the channels
		for(int i=0; i<w; i++) {
			for (int k=0;k<bpp;k++) {
				ogf_swap(row1[bpp*i+k], row2[bpp*i+k]);
			}
		}
	}
	return true;
}

bool ImageSerializer::swap_channels(Image& image, 
									int channel0, int channel1) const 
{
	int bpp=image.bytes_per_pixel();
	int h = image.height() ;
	int w = image.width() ;
	int row_len = w * bpp ;

	if (bpp<=channel0 || bpp<=channel1) {
		return false;
	}

	for(int j=0; j< h; j++) {
		// get a pointer to the line we will swap channels on
		Memory::pointer row1 = image.base_mem() + j * row_len ;
		// for each point on line, swap the two channels
		for(int i=0; i<w; i++) {
			ogf_swap(row1[bpp*i+channel0], row1[bpp*i+channel1]);
		}
	}
	return true;
}

bool ImageSerializer::rgb_to_bgr(Image& image) const{
	if (image.bytes_per_pixel()<3) {
		return false;
	}
	return swap_channels(image,0,2);
}

