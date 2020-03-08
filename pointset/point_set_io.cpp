#include "point_set_io.h"
#include "point_set_serializer_ply.h"
#include "point_set.h"
#include "../basic/stop_watch.h"
#include "../basic/file_utils.h"
#include "../basic/logger.h"

#include <fstream>
#include <list>


PointSet* PointSetIO::read(const std::string& file_name)
{
	std::ifstream in(file_name.c_str()) ;
	if(in.fail()) {
		Logger::err(title()) << "cannot open file: " << file_name << std::endl;
		return nil;
	}
	in.close();

	Logger::out(title()) << "reading file..." << std::endl;

	StopWatch w;
	PointSet* object = nil;

	std::string ext = FileUtils::extension(file_name);
	String::to_lowercase(ext);

	if (ext == "ply") {
		object = PointSetSerializer_ply::load(file_name);
	} 

	else {
		PointSet* point_set = new PointSet;
		if (ext == "pnc") 
			load_pnc(point_set, file_name);
		else if (ext == "bpnc") 
			load_bpnc(point_set, file_name);

		else {
			Logger::err(title()) << "reading file failed (unknown file format)" << std::endl;
			delete point_set;
			return nil;
		}
		
		if (point_set->num_points() < 1) {
			Logger::err(title()) << "reading file failed (no data exist)" << std::endl;
			delete point_set;
			return nil;
		}
		object = point_set;
	}	

	Logger::out(title()) << "reading file done. Time: "
		<< w.elapsed() << " seconds" << std::endl;
	
	return object;
}

bool PointSetIO::save(const std::string& file_name, const PointSet* point_set) {
	if (!point_set) {
		Logger::err(title()) << "Point set is null" << std::endl;
		return false;
	}
	
	std::ofstream out(file_name.c_str()) ;
	if(out.fail()) {
		Logger::err(title()) << "cannot open file: \'" << file_name << "\' for writing" << std::endl;
		return false;
	}
	Logger::out(title()) << "saving file..." << std::endl;
	out.close();

	StopWatch w;
	std::string ext = FileUtils::extension(file_name);
	String::to_lowercase(ext);
	
	out.precision(16);

 	if (ext == "ply")
 		PointSetSerializer_ply::save(file_name, point_set);
	else if (ext == "pnc")
		save_pnc(point_set, file_name);
	else if (ext == "bpnc")
		save_bpnc(point_set, file_name);
	else {
		Logger::err(title()) << "saving file failed (unknown file format)" << std::endl;
		return false;
	}

	Logger::out(title()) << "saving file done. Time: "
		<< w.elapsed() << " seconds" << std::endl;

	return true;
}


void PointSetIO::load_pnc(PointSet* pointSet, const std::string& file_name) {
	std::ifstream input(file_name.c_str()) ;
	if(input.fail()) {
		Logger::err(title()) << "could not open file\'" << file_name << "\'" << std::endl;
		return ;
	}

	std::vector<vec3f>& points = pointSet->points();	points.clear();
	std::vector<vec3f>& normals = pointSet->normals();	normals.clear();
	std::vector<vec3f>& colors = pointSet->colors();	colors.clear();

	float x, y, z, nx, ny, nz, r, g, b;
	while(!input.eof()) {
		input >> x >> y >> z >> nx >> ny >> nz >> r >> g >> b;

		if (!input.fail()) {
			points.push_back(vec3f(x, y, z));
			normals.push_back(vec3f(nx, ny, nz));
			colors.push_back(vec3f(r, g, b));
		}
	}
}


void PointSetIO::save_pnc(const PointSet* pointSet, const std::string& file_name) {
	// open file
	std::ofstream output(file_name.c_str()) ;
	if(output.fail()) {
		Logger::err(title()) << "could not open file\'" << file_name << "\'" << std::endl;
		return ;
	}
	output.precision(16);

	int num = pointSet->num_points();
	float* points = &(const_cast<PointSet*>(pointSet)->points()[0].x);
	float* normals = 0;
	float* colors = 0;

	if (pointSet->has_normals())
		normals = &(const_cast<PointSet*>(pointSet)->normals()[0].x);
	else {
		normals = new float[num * 3];
		memset(normals, 0, num * 3 * 4);
	}

	if (pointSet->has_colors())
		colors = &(const_cast<PointSet*>(pointSet)->colors()[0].x);
	else {
		colors = new float[num * 3];
		memset(colors, 0, num * 3 * 4);
	}

	for (int i = 0; i < num; ++i) {
		output 
			<< points[i * 3] << " " << points[i * 3 + 1] << " " << points[i * 3 + 2] << " "
			<< normals[i * 3] << " " << normals[i * 3 + 1] << " " << normals[i * 3 + 2] << " "
			<< colors[i * 3] << " " << colors[i * 3 + 1] << " " << colors[i * 3 + 2] << std::endl;
	}
}


// each line: x y z nx ny nz r g b. All are float
void PointSetIO::load_bpnc(PointSet* pointSet, const std::string& file_name) {
	std::ifstream input(file_name.c_str(), std::fstream::binary) ;
	if(input.fail()) {
		Logger::err(title()) << "could not open file\'" << file_name << "\'" << std::endl;
		return ;
	}

	// check size of types
	int line_size = sizeof(float) * 9;

	std::streamoff begin_pos = input.tellg();
	input.seekg(0, std::ios::end);
	std::streamoff end_pos = input.tellg();
	// num of points in the file
	int num = static_cast<int>(end_pos - begin_pos) / line_size;  

	input.seekg(0, std::ios::beg); 

	float* data = new float[num * 9]; 
	input.read((char*)data, num * line_size);	// read the entire blocks

	std::vector<vec3f>& points = pointSet->points();	points.resize(num);
	std::vector<vec3f>& normals = pointSet->normals();	normals.resize(num);
	std::vector<vec3f>& colors = pointSet->colors();	colors.resize(num);

	for (int i = 0; i<num; ++i) {
		float* p = data + i*9; 
		points[i] = vec3f(p);

		float* n = data + i * 9 + 3;
		normals[i] = vec3f(n);

		float* c = data + i * 9 + 6;
		colors[i] = vec3f(c);
	}
	delete [] data;
}


void PointSetIO::save_bpnc(const PointSet* pointSet, const std::string& file_name) {
	int line_size = sizeof(float) * 9;

	// open file
	std::ofstream output(file_name.c_str(), std::fstream::binary) ;
	if(output.fail()) {
		Logger::err(title()) << "could not open file\'" << file_name << "\'" << std::endl;
		return ;
	}

	int num = pointSet->num_points();
	float* points = &(const_cast<PointSet*>(pointSet)->points()[0].x);
	float* normals = 0;
	float* colors = 0;

	if (pointSet->has_normals())
		normals = &(const_cast<PointSet*>(pointSet)->normals()[0].x);
	else {
		normals = new float[num * 3];
		memset(normals, 0, num * 3 * 4);
	}

	if (pointSet->has_colors())
		colors = &(const_cast<PointSet*>(pointSet)->colors()[0].x);
	else {
		colors = new float[num * 3];
		memset(colors, 0, num * 3 * 4);
	}

	for (int i = 0; i < num; ++i) {
		float* pt = points + (i * 3);	output.write((char*)pt, 12);
		float* nm = normals + (i * 3);	output.write((char*)nm, 12);
		float* cl = colors + (i * 3);	output.write((char*)cl, 12);
	}

	if (!pointSet->has_normals())
		delete[] normals;
	if (!pointSet->has_colors())
		delete[] colors;
}
