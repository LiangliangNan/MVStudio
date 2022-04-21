
#ifndef __IMAGE_IO__
#define __IMAGE_IO__



#include <string>


class Image ;
class ImageSerializer ;

class ImageIO
{
public:
	static Image*	read(const std::string& file_name);
	static bool		save(const std::string& file_name, const Image* image) ;

	// only read the image file header
	static bool		query_image_size(const std::string& file_name, int& width, int& height);

	static ImageSerializer* resolve_serializer(const std::string& file_name) ;

};


#endif

