#include "image_data.h"
#include "../basic/basic_types.h"
#include "../image/image.h"
#include "../image/image_io.h"


namespace sfm {

	ImageData::ImageData()
		: image_(nil)
		, ignore_in_bundle(false)
		, image_loaded(false)
		, keys_loaded(false)
		, width_(0)
		, height_(0)
		, cached_size_(false)
		, num_keys_(0)
		, cached_keys_(false)
		, has_init_focal(false)
		, init_focal(0.0)
		, known_intrinsics(false)
	{
	}

	ImageData::~ImageData() {
		if (image_)
			unload_image();
	}

	int ImageData::width() {
		if (image_loaded)
			return (image_->width());
		else {
			if (cached_size_)
				return width_;

			int w = -1;
			int h = -1;
			if (ImageIO::query_image_size(image_file, w, h)) {
				width_ = w;
				height_ = h;
				cached_size_ = true;
				return width_;
			}
			return -1; // error 
		}
	}


	int ImageData::height() {
		if (image_loaded)
			return (image_->height());
		else {
			if (cached_size_)
				return height_;

			int w = -1;
			int h = -1;
			if (ImageIO::query_image_size(image_file, w, h)) {
				width_ = w;
				height_ = h;
				cached_size_ = true;
				return height_;
			}
			return -1; // error 
		}
	}

	void ImageData::load_image() {
		if (image_ && image_loaded)
			return;

		image_ = ImageIO::read(image_file);
		if (image_)
			image_loaded = true;
	}

	void ImageData::unload_image() {
		if (!image_loaded || !image_)
			return;
		
		delete image_;
		image_ = nil;

		image_loaded = false;
	}

	int  ImageData::num_of_keys() {
		if (keys_loaded) {
			return (int)keys.size();
		}
		else {
			if (!cached_keys_)
				num_keys_ = read_num_of_keys(key_file);

			return num_keys_;
		}
	}

	void ImageData::load_keys() {
		if (keys_loaded)
			return;   /* Already loaded the keys */

		// Try to read data from a key point file
		std::vector<Keypoint> kps;
		//if (!read_keys_ascii(key_file, kps))
		if (!read_keys_binary(key_file, kps))
				return;

		// Flip y-axis to make things easier
		for (int i = 0; i < (int)kps.size(); i++) {
			kps[i].y = height() - kps[i].y - 1.0f;
		}

		// Now make the image center the origin
		for (int i = 0; i < (int)kps.size(); i++) {
			kps[i].x -= 0.5f * (width() - 1);
			kps[i].y -= 0.5f * (height() - 1);
		}

		keys = kps;
		keys_loaded = true;
	}

	void ImageData::unload_keys() {
		keys.clear();
		keys_loaded = false;
	}

	void ImageData::load_key_colors() {
		bool unload = false;
		if (!image_loaded) {
			load_image();
			unload = true;
		}

		int w = width();
		int h = height();

		int num_keys = (int)keys.size();

		for (int i = 0; i < num_keys; i++) {
			double x = keys[i].x;
			double y = keys[i].y;

			Colorf c = pixel_interpolate(image_, x + 0.5 * w, y + 0.5 * h);

			keys[i].r = iround(c.r());
			keys[i].g = iround(c.g());
			keys[i].b = iround(c.b());
		}

		if (unload)
			unload_image();
	}
}