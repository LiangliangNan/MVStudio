
#ifndef ___IMAGE_SERIALIZER_BMP__
#define ___IMAGE_SERIALIZER_BMP__

#include "image_common.h"
#include "image_serializer.h"

class IMAGE_API ImageSerializer_bmp : public ImageSerializer {
public:
	virtual Image*	serialize_read(std::istream& in) ;
	virtual bool	serialize_write(std::ostream& out,const Image *image);

	virtual bool	read_supported() const {return true ;}
	virtual bool	write_supported() const {return true ;}

	virtual bool	binary() const {return true;};
	virtual bool	streams_supported() const {return true ;}

	virtual	bool	query_image_size(const std::string& file_name, int& width, int& height);
} ;


#endif

