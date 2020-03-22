
#ifndef _DENSE_RECONSTRUCTION_H_
#define _DENSE_RECONSTRUCTION_H_


#include "sparse_reconstruction.h"


#include <string>
#include <vector>


class PointSet;
class Project;

class DenseReconstruction {
public:
	DenseReconstruction(Project* proj);
	virtual ~DenseReconstruction();

	static std::string title() { return "DenseRecon"; }

	bool run_pmvs(PointSet* pset = nil);

private:
	bool convert_bundler_to_pmvs();

private:
	Project* project_;
} ;



#endif

