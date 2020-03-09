#ifndef __IMAGE_SERIALIZER_JPEG__
#define __IMAGE_SERIALIZER_JPEG__


#include "image_serializer.h"


class ImageSerializer_jpeg : public ImageSerializer {
public:
	virtual Image*	serialize_read(const std::string& file_name) ;
	virtual bool	serialize_write(
		const std::string& file_name, const Image* image
		);

	virtual bool	read_supported() const ;
	virtual bool	streams_supported() const ;

	virtual	bool	query_image_size(const std::string& file_name, int& width, int& height);
} ;

#endif

