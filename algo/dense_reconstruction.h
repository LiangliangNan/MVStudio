
#ifndef _DENSE_RECONSTRUCTION_H_
#define _DENSE_RECONSTRUCTION_H_

#include "sparse_reconstruction.h"


#include <string>
#include <vector>

namespace easy3d {
    class PointCloud;
}

class Project;

class  DenseReconstruction {
public:
	DenseReconstruction(Project* proj);
	virtual ~DenseReconstruction();

	bool run_pmvs(easy3d::PointCloud* pset = nullptr);

private:
	bool convert_bundler_to_pmvs();
	bool convert_pba_to_pmvs();

private:
	Project* project_;
} ;



#endif

