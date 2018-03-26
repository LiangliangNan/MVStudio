#include "keys.h"

#include <fstream>
#include <iostream>



namespace sfm {

	/* Returns the number of keys in a file */
	int read_num_of_keys(const std::string& filename) {
		std::ifstream input(filename.c_str(), std::fstream::binary);
		if (input.fail()) { // so could be ascii file
			input.open(filename.c_str());
			if (input.fail()) {
				std::cerr << "could not open file: \'" << filename << "\'" << std::endl;
				return 0;
			}	
			int num, len;
			input >> num >> len;
			if (!input.fail()) {
				return num;
			}
			else {
				std::cerr << "Invalid key point file." << std::endl;
				return 0;
			}
		}
		else
		{
			int num;
			input.read((char*)(&num), sizeof(int));
			if (!input.fail()) {
				return num;
			}
			else {
				std::cerr << "Invalid key point file." << std::endl;
				return 0;
			}
		}
	}

	/* Read keypoints from the given file pointer and return the list of
	* keypoints.  The file format starts with 2 integers giving the total
	* number of keypoints and the size of descriptor vector for each
	* keypoint (currently assumed to be 128). Then each keypoint is
	* specified by 4 floating point numbers giving subpixel row and
	* column location, scale, and orientation (in radians from -PI to
	* PI).  Then the descriptor vector for each keypoint is given as a
	* list of integers in range [0,255]. */
	bool read_keys_ascii(const std::string& filename, std::vector<Keypoint>& keys) {
		std::ifstream input(filename.c_str());
		if (input.fail()) {
			std::cerr << "could not open file: \'" << filename << "\'" << std::endl;
			return false;
		}
		int num, len;
		input >> num >> len;
		if (input.fail()) {
			std::cerr << "Invalid key point file." << std::endl;
			return false;
		}

		if (len != 128) {
			std::cerr << "Keypoint descriptor length invalid (should be 128)" << std::endl;
			return false;
		}

		keys.resize(num);
		for (int i = 0; i < num; i++) {
			float x, y, scale, orient;
			input >> y >> x >> scale >> orient;	// row first
			if (input.fail()) {
				std::cerr << "Invalid key point file format" << std::endl;
				return false;
			}

			int desc[128];
			for (int j = 0; j < len; j++) {
				int val;
				input >> val;
				if (input.fail() || val < 0 || val > 255) {
					std::cerr << "Invalid key point file value" << std::endl;
					return false;
				}
				desc[j] = val;
			}

			keys[i].x = x;
			keys[i].y = y;

			// add code here to take also the 'scale, orient, desc' ...
			// ...
		}

		return true;
	}

	/* Read keys from binary file */
	bool read_keys_binary(const std::string& filename, std::vector<Keypoint>& keys) {
		std::ifstream input(filename.c_str(), std::fstream::binary);
		if (input.fail()) {
			std::cerr << "could not open file: \'" << filename << "\'" << std::endl;
			return false;
		}

		int num, dim;
		input.read((char*)(&num), sizeof(int));
		input.read((char*)(&dim), sizeof(int));

		keys.resize(num);

		typedef struct {
			float y;  // row first
			float x;
			float scale;
			float orient;
			int	desc[128];
		} keypt_t;
		keypt_t* info = new keypt_t[num];

		input.read((char*)info, sizeof(keypt_t) * num);
		for (int i = 0; i < num; i++) {
			keys[i].x = info[i].x;
			keys[i].y = info[i].y;

			// add code here to take also the 'scale, orient, desc' ...
			// ...
		}

		delete[] info;

		return true;
	}

}