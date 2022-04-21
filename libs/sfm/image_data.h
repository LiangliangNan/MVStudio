#ifndef _SFM_IMAGE_DATA_H_
#define _SFM_IMAGE_DATA_H_

#include "keys.h"

#include "camera.h"

#include <string>
#include <vector>


class Image;

namespace sfm {

	class ImageData
	{
	public:
		ImageData();
		~ImageData();

		// Return the image size 
		int  width();
		int  height();

		void load_image();
		void unload_image();

		int  num_of_keys();

		void load_keys();
		void unload_keys();

		// Read key colors
		void load_key_colors();

	public:
		std::vector<Keypoint>	keys;			// Key points in this image
		std::vector<bool>		key_flags;
		std::vector<int>		visible_points;	// Indices of points visible in this image
		std::vector<int>		visible_keys;
		Camera			camera;		// Information on the camera used to capture this image
		
		/* Radial distortion parameters */
		bool			known_intrinsics;
		double		K[9];
		double		k[5];

		bool			ignore_in_bundle;  /* Ignore this image during bundle adjustment */

		bool			image_loaded;
		bool			keys_loaded;
		std::string		image_file;
		std::string		key_file;

		Image*		image_;

		// Focal length parameters 
		bool			has_init_focal;	// Do we have an initial focal length? 
		double		init_focal;		// Initial focal length 

	private:
		int			width_;
		int			height_;		// Cached dimensions
		bool			cached_size_;

		int			num_keys_;      // Cached number of keys
		bool			cached_keys_;

	};

}


#endif /* __IMAGE_DATA_H__ */
