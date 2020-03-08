
#ifndef __IMAGE_SERIALIZER_PNG__
#define __IMAGE_SERIALIZER_PNG__

#include "image_common.h"
#include "image_serializer.h"


class IMAGE_API ImageSerializer_png : public ImageSerializer {
public:
	virtual Image*	serialize_read(const std::string& file_name) ;
	virtual bool	read_supported() const ;
	virtual bool	serialize_write(
		const std::string& file_name, const Image* image
		) ;

	virtual bool	write_supported() const ;
	virtual bool	streams_supported() const ;
	virtual bool	binary() const ;

	virtual	bool	query_image_size(const std::string& file_name, int& width, int& height);
} ;


#endif

