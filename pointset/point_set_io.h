#ifndef _POINT_SET_IO_H_
#define _POINT_SET_IO_H_



#include <string>


class PointSet;

class PointSetIO
{
public:
	static std::string title() { return "PointSetIO"; }

	// for both point cloud and mesh
	static PointSet* read(const std::string& file_name);

	// save the point set to a file. return false if failed.
	static bool		 save(const std::string& file_name, const PointSet* point_set);

protected:
	// each line with point, normal and color: (x, y, z, nx, ny, nz, r, g, b)
	static void load_pnc(PointSet* pointSet, const std::string& file_name);
	static void save_pnc(const PointSet* pointSet, const std::string& file_name);
	static void load_bpnc(PointSet* pointSet, const std::string& file_name);
	static void save_bpnc(const PointSet* pointSet, const std::string& file_name);
};

#endif