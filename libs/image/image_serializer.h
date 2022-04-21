
#ifndef _IMAGE_IO_IMAGE_SERIALIZER_H_
#define _IMAGE_IO_IMAGE_SERIALIZER_H_


#include "../basic/counted.h"

#include <iostream>


class Image ;

class ImageSerializer : public Counted {
public:

	virtual Image*	serialize_read(const std::string& file_name) ;
	virtual bool	serialize_write(
		const std::string& file_name, const Image* image
		) ;

	virtual Image*	serialize_read(std::istream& stream) ;
	virtual bool	serialize_write(
		std::ostream& stream, const Image* image
		) ;

	virtual	bool	query_image_size(const std::string& file_name, int& width, int& height) = 0;

	/**
	* checks whether the stream should be opened
	* in text or binary mode. Default returns true.
	*/
	virtual bool binary() const ;

	/**
	* checks whether reading and writing to streams is
	* supported.
	*/
	virtual bool streams_supported() const ;

	/**
	* checks whether reading is implemented.
	*/
	virtual bool read_supported() const ;

	/**
	* checks whether writing is implemented.
	*/
	virtual bool write_supported() const ;

	/**
	* function provided to flip an image vertically
	*/
	virtual bool flip_image(Image& image) const;

	/**
	* function provided to swap channels on image
	*/
	virtual bool swap_channels(Image& image, int channel0, int channel1) const;

	/**
	* function provided to convert rgb image to bgr images 
	* (with or without alpha channel)
	*/
	virtual bool rgb_to_bgr(Image& image) const;
} ; 

typedef SmartPointer<ImageSerializer> ImageSerializer_var ;


#endif

