
#ifndef _SFM_KEYS_H_
#define _SFM_KEYS_H_

#include "../basic/basic_types.h"

#include <vector>
#include <string>


namespace sfm {

	class Keypoint {
	public:
		Keypoint() : x(0), y(0), r(0), g(0), b(0), extra(-1), track(-1) {}
		Keypoint(float _x, float _y) : x(_x), y(_y), r(0), g(0), b(0), extra(-1), track(-1) {}

		virtual ~Keypoint() {}

		virtual unsigned char *desc() { return nil; }

		float x;
		float y;			/* Subpixel location of keypoint. */

		unsigned char r;
		unsigned char g;
		unsigned char b;	/* Color of this key */

		int extra;  /* 4 bytes of extra storage */
		int track;  /* Track index this point corresponds to */
	};



	/* Data struct for matches */
	class KeypointMatch {
	public:
		KeypointMatch() { }
		KeypointMatch(int idx1, int idx2) : key_idx1(idx1), key_idx2(idx2) { }

		int key_idx1, key_idx2;
	};


	/* Returns the number of keys in a file */
	int read_num_of_keys(const std::string& filename);

	/* Read keypoints from the given file name and return the list of
	 * keypoints.  The file format starts with 2 integers giving the total
	 * number of keypoints and the size of descriptor vector for each
	 * keypoint (currently assumed to be 128). Then each keypoint is
	 * specified by 4 floating point numbers giving subpixel row and
	 * column location, scale, and orientation (in radians from -PI to
	 * PI).  Then the descriptor vector for each keypoint is given as a
	 * list of integers in range [0,255]. */
	bool read_keys_ascii(const std::string& filename, std::vector<Keypoint>& keys);
	bool read_keys_binary(const std::string& filename, std::vector<Keypoint>& keys);
}


#endif 
