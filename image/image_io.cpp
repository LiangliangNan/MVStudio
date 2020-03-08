
#include "image_io.h"
#include "image.h"
#include "image_serializer.h"
#include "image_serializer_png.h"
#include "image_serializer_jpeg.h"
#include <easy3d/util/logging.h>
#include <easy3d/util/file_system.h>


using namespace easy3d;


Image* ImageIO::read(const std::string& file_name) {
	ImageSerializer* serializer = resolve_serializer(file_name);
	if (!serializer) {
        Image* img = serializer->serialize_read(file_name);
        delete serializer;
        return img;
    }
	return nullptr;
}


bool ImageIO::save(const std::string& file_name, const Image* image) {
	ImageSerializer* serializer = resolve_serializer(file_name);
	if (!serializer) {
        bool result = serializer->serialize_write(file_name, image);
        delete serializer;
        return result;
    }

	return false;
}

ImageSerializer* ImageIO::resolve_serializer(const std::string& file_name) {
	std::string extension = file_system::extension(file_name) ;
	if(extension.length() == 0) {
		LOG(ERROR) << "No extension in file name" << std::endl ;
		return nullptr ;
	}

	ImageSerializer* serializer = nullptr;

	if ( extension == "png" )
		serializer = new ImageSerializer_png();
	else if ( extension == "jpg" )
		serializer = new ImageSerializer_jpeg();
	else { 	
		LOG(ERROR) << "Unknown image file format \'" << extension << "\'" << std::endl;
		return nullptr;
	}

	return serializer;
}


bool ImageIO::query_image_size(const std::string& file_name, int& width, int& height) {
	ImageSerializer* serializer = resolve_serializer(file_name);
	if (!serializer) {
	    bool result = serializer->query_image_size(file_name, width, height);
	    delete serializer;
        return result;
	}

	return false;
}
