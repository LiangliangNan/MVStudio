
#ifndef _IMAGE_MATCHING_H_
#define _IMAGE_MATCHING_H_




#include <string>
#include <vector>


class Project;

class ImageMatching {
public:
	ImageMatching(Project* proj);
	virtual ~ImageMatching();

	static std::string title() { return "ImageMatching"; }

	void apply();

private:
	void extract_key_points();
	void match_key_points();

private:
	Project*	project_;
} ;



#endif

