
#include "point_set_serializer_ply.h"
#include "rply.h"
#include "point_set.h"
#include "../basic/logger.h"
#include "../basic/basic_types.h"

#include <cassert>



class PlyPointSetLoad 
{
public:
	PlyPointSetLoad(PointSet* pset)
		: points(pset->points())
		, normals(pset->normals())
		, colors(pset->colors())
	{
		points.clear();
		normals.clear();
		colors.clear();
	}

	bool load(const std::string& filename) {
		p_ply ply = ply_open(filename.c_str(), nil, 0, nil) ;
		if(ply == nil) {
			Logger::err("PlyPointSetLoad") << filename << ": could not open" << std::endl;
			return false ;
		}

		if(!ply_read_header(ply)) {
			Logger::err("PlyPointSetLoad") << filename << ": invalid PLY file" << std::endl;
			ply_close(ply) ;
			return false ;
		}

		check_for_colors_and_normals(ply) ;

		long nvertices = ply_set_read_cb(ply, "vertex", "x", PlyPointSetLoad::vertex_cb, this, 0) ;
		ply_set_read_cb(ply, "vertex", "y", PlyPointSetLoad::vertex_cb, this, 1) ;
		ply_set_read_cb(ply, "vertex", "z", PlyPointSetLoad::vertex_cb, this, 2) ;

		if(!ply_read(ply)) {
			Logger::err("PlyPointSetLoad") 
				<< filename << ": problem occurred while parsing PLY file" << std::endl;
		}

		ply_close(ply) ;
		return end_point_set() ;
	}

protected:
	void check_for_colors_and_normals(p_ply ply) {
		p_ply_element element = nil;

		std::string str_red, str_green, str_blue;

		bool has_normals = false;
		bool has_normals_Neil = false;

		for (;;) {
			element = ply_get_next_element(ply, element);
			if (element == nil) { break; }
			const char* elt_name = nil;
			ply_get_element_info(element, &elt_name, nil);

			if (!strcmp(elt_name, "vertex")) {
				p_ply_property property = nil;
				for (;;) {
					property = ply_get_next_property(element, property);
					if (property == nil)
						break;

					const char* prop_name = nil;
					ply_get_property_info(property, &prop_name, nil, nil, nil);
					if (!strcmp(prop_name, "r"))	str_red = "r";
					if (!strcmp(prop_name, "g"))	str_green = "g";
					if (!strcmp(prop_name, "b"))	str_blue = "b";
					if (!strcmp(prop_name, "red"))		str_red = "red";
					if (!strcmp(prop_name, "green"))	str_green = "green";
					if (!strcmp(prop_name, "blue"))		str_blue = "blue";
					if (!strcmp(prop_name, "diffuse_red"))		str_red = "diffuse_red";
					if (!strcmp(prop_name, "diffuse_green"))	str_green = "diffuse_green";
					if (!strcmp(prop_name, "diffuse_blue"))		str_blue = "diffuse_blue";

					has_normals = has_normals || !strcmp(prop_name, "nx");
					has_normals = has_normals || !strcmp(prop_name, "ny");
					has_normals = has_normals || !strcmp(prop_name, "nz");
					has_normals_Neil = has_normals || !strcmp(prop_name, "vsfm_cnx");  // for Neil Smith's
					has_normals_Neil = has_normals || !strcmp(prop_name, "vsfm_cny");
					has_normals_Neil = has_normals || !strcmp(prop_name, "vsfm_cnz");
				}
			}
		}

		if (str_red == "r" && str_green == "g" && str_blue == "b") {
			has_colors_ = true;
			color_mult_ = 1.0f;
		}
		else if (
			(str_red == "red" && str_green == "green" && str_blue == "blue") ||
			(str_red == "diffuse_red" && str_green == "diffuse_green" && str_blue == "diffuse_blue")
			)
		{
			has_colors_ = true;
			color_mult_ = 1.0f / 255.0f;
		}
		else {
			has_colors_ = false;
		}

		if (has_colors_) {
			ply_set_read_cb(ply, "vertex", str_red.c_str(), PlyPointSetLoad::color_cb, this, 0);
			ply_set_read_cb(ply, "vertex", str_green.c_str(), PlyPointSetLoad::color_cb, this, 1);
			ply_set_read_cb(ply, "vertex", str_blue.c_str(), PlyPointSetLoad::color_cb, this, 2);
		}

		has_normals_ = has_normals;
		if (has_normals) {
			has_normals_ = true;
			ply_set_read_cb(ply, "vertex", "nx", PlyPointSetLoad::normal_cb, this, 0);
			ply_set_read_cb(ply, "vertex", "ny", PlyPointSetLoad::normal_cb, this, 1);
			ply_set_read_cb(ply, "vertex", "nz", PlyPointSetLoad::normal_cb, this, 2);
		}
		else
			has_normals_ = false;
	}

	static PlyPointSetLoad* plyload(p_ply_argument argument) {
		PlyPointSetLoad* result = nil ;
		ply_get_argument_user_data(argument, (void**)(&result), nil) ;
		ogf_assert(result != nil) ;
		return result ;
	}

	static int vertex_cb(p_ply_argument argument) {
		return plyload(argument)->add_vertex_data(argument) ;
	}

	static int color_cb(p_ply_argument argument) {
		return plyload(argument)->add_color_data(argument) ;
	}

	static int normal_cb(p_ply_argument argument) {
		return plyload(argument)->add_normal_data(argument) ;
	}

	int add_vertex_data(p_ply_argument argument) {
		long coord ;
		ply_get_argument_user_data(argument, nil, &coord);
		ogf_assert(coord >= 0 && coord < 3) ;
		xyz_[coord] = float(ply_get_argument_value(argument)) ;
		if(coord == 2) { 
			points.push_back(vec3f(xyz_));
		}
		return 1;
	}

	int add_color_data(p_ply_argument argument) {
		long coord ;
		ply_get_argument_user_data(argument, nil, &coord);
		ogf_assert(coord >= 0 && coord < 3) ;
		rgb_[coord] = float(ply_get_argument_value(argument)) * color_mult_ ;
		if(coord == 2) { 
			colors.push_back(vec3f(rgb_));
		}
		return 1 ;
	}

	int add_normal_data(p_ply_argument argument) {
		long coord ;
		ply_get_argument_user_data(argument, nil, &coord);
		ogf_assert(coord >= 0 && coord < 3) ;
		normal_[coord] = float(ply_get_argument_value(argument));
		if(coord == 2) { 
			normals.push_back(vec3f(normal_));
		}
		return 1 ;
	}

	//////////////////////////////////////////////////////////////////////////

	bool end_point_set() {
		return points.size() > 0;
	}

protected:
	std::vector<vec3f>& points;
	std::vector<vec3f>& normals;
	std::vector<vec3f>& colors;

	float			xyz_[3];
	float			rgb_[3];
	float			normal_[3];

	bool			has_normals_;
	bool			has_colors_ ;
	float			color_mult_;
} ;

//__________________________________________________________

PointSet*	PointSetSerializer_ply::load(const std::string& file_name) {
	PointSet* pset = new PointSet;
	PlyPointSetLoad loader(pset) ;
	if (loader.load(file_name))
		return pset;
	else {
		delete pset;
		return nil;
	}
}


class PlyPointSetSave {
public:
	PlyPointSetSave(const PointSet* obj) 
		: object_(obj)
		, color_mult_(255.0)
	{ }

	bool save(const std::string& filename) {
		p_ply ply = ply_create(filename.c_str(), PLY_LITTLE_ENDIAN, nil, 0, nil) ;

		if(ply == nil) {
			Logger::err("PlyPointSetSave") << filename << ": could not open" << std::endl;
			return false ;
		}

		//////////////////////////////////////////////////////////////////////////

		if (!ply_add_comment(ply, "saved by liangliang.nan@gmail.com")) {
			Logger::err("PlyPointSetSave") << "unable to add comment" << std::endl;
			ply_close(ply) ;
			return false ;
		}

		int num = object_->num_points();
		if (!ply_add_element(ply, "vertex", num)) {
			Logger::err("PlyPointSetSave") << "unable to add element \'vertex\'" << std::endl;
			ply_close(ply) ;
			return false ;
		}

		e_ply_type length_type, value_type;
		length_type = value_type = static_cast<e_ply_type>(-1);
		std::string pos[3] = { "x", "y", "z" };
		for (unsigned int i=0; i<3; ++i) {
			if (!ply_add_property(ply, pos[i].c_str(), PLY_FLOAT, length_type, value_type)) {
				Logger::err("PlyPointSetSave") << "unable to add property \'" << pos[i] << "\'" << std::endl;
				ply_close(ply) ;
				return false ;
			}
		}
	
		if (object_->has_normals()) {
			std::string normal[3] = { "nx", "ny", "nz" };
			for (unsigned int i = 0; i < 3; ++i) {
				if (!ply_add_property(ply, normal[i].c_str(), PLY_FLOAT, length_type, value_type)) {
					Logger::err("PlyPointSetSave") << "unable to add property \'" << pos[i] << "\'" << std::endl;
					ply_close(ply);
					return false;
				}
			}
		}

		if (object_->has_colors()) {
			std::string color[4] = { "diffuse_red", "diffuse_green", "diffuse_blue" };
			for (unsigned int i = 0; i < 3; ++i) {
				if (!ply_add_property(ply, color[i].c_str(), PLY_UCHAR, length_type, value_type)) {
					Logger::err("PlyPointSetSave") << "unable to add property \'" << color[i] << "\'" << std::endl;
					ply_close(ply);
					return false;
				}
			}
		}

		if(!ply_write_header(ply)) {
			Logger::err("PlyPointSetSave") << filename << ": invalid PLY file" << std::endl;
			ply_close(ply) ;
			return false ;
		}

		//////////////////////////////////////////////////////////////////////////

		const std::vector<vec3f>& points = object_->points();
		const std::vector<vec3f>& normals = object_->normals();
		const std::vector<vec3f>& colors = object_->colors();

		for (int i = 0; i < num; ++i) {
			ply_write(ply, points[i].x);
			ply_write(ply, points[i].y);
			ply_write(ply, points[i].z);

			if (object_->has_normals()) {
				ply_write(ply, normals[i].x);
				ply_write(ply, normals[i].y);
				ply_write(ply, normals[i].z);
			}

			if (object_->has_colors()) {
				float r = colors[i].x * color_mult_;	ogf_clamp(r, 0.0f, 255.0f);
				float g = colors[i].y * color_mult_;	ogf_clamp(g, 0.0f, 255.0f);
				float b = colors[i].z * color_mult_;	ogf_clamp(b, 0.0f, 255.0f);
				ply_write(ply, r);
				ply_write(ply, g);
				ply_write(ply, b);
			}
		}
		
		ply_close(ply);
		return true ;
	}

protected:
	const PointSet*	object_;
	float			color_mult_;
} ;


bool PointSetSerializer_ply::save(const std::string& file_name, const PointSet* obj) {
	PlyPointSetSave plysave(obj) ;
	return  plysave.save(file_name) ;
}