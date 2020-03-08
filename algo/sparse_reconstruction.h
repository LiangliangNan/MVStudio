
#ifndef _SPARSE_RECONSTRUCTION_H_
#define _SPARSE_RECONSTRUCTION_H_


#include "../sfm/sfm.h"


#include <string>
#include <vector>


namespace easy3d {
    class PointCloud;
}

class Project;

class  SparseReconstruction {
public:
	SparseReconstruction(Project* proj);
	virtual ~SparseReconstruction();

	static std::string title() { return "SparseRecon"; }

	bool apply(easy3d::PointCloud* pset = nullptr);

private:
	//------------------ SfM ------------------
	void run_sfm(easy3d::PointCloud* pset);

	//------------------ MVS ------------------
	void convert_bundler_to_pmvs();
	void convert_pba_to_pmvs();
	void run_mvs();

private:
	Project*	project_;
} ;



#endif

