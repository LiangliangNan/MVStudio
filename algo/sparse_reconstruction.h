
#ifndef _SPARSE_RECONSTRUCTION_H_
#define _SPARSE_RECONSTRUCTION_H_

#include "algo_common.h"

#include "../sfm/sfm.h"


#include <string>
#include <vector>


class PointSet;
class Project;

class ALGO_API SparseReconstruction {
public:
	SparseReconstruction(Project* proj);
	virtual ~SparseReconstruction();

	static std::string title() { return "SparseRecon"; }

	bool apply(PointSet* pset = nil);

private:
	//------------------ SfM ------------------
	void run_sfm(PointSet* pset);

	//------------------ MVS ------------------
	void convert_bundler_to_pmvs();
	void convert_pba_to_pmvs();
	void run_mvs();

private:
	Project*	project_;
} ;



#endif

