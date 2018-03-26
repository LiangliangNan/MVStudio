#include "sfm.h"
#include "../basic/logger.h"
#include "../basic/file_utils.h"
#include "../math/matrix_driver.h"

#include <fstream>


namespace sfm {

	/* Load a list of image names from a file */
	void SfM::load_image_names() {
		std::ifstream input(option_.list_file.c_str());
		if (input.fail()) {
			Logger::warn(title()) << "could not open list file \'" << option_.list_file << "\'" << std::endl;
			return;
		}

		image_data_.clear();

		while (!input.eof()) {
			std::string line;
			getline(input, line);
			if (input.fail())
				break;

			std::istringstream stream(line);
			if (stream.fail())
				break;

			std::string image_file, dummy;
			int img_width, img_height;
			float focal_length;
			stream >> image_file >> dummy >> img_width >> dummy >> img_height >> dummy >> focal_length;

			if (stream.fail())
				break;

			if (!FileUtils::is_file(image_file)) {
				Logger::warn(title()) << "image file \'" << FileUtils::simple_name(image_file) << "\' doesn't exist" << std::endl;
				continue;
			}
				
			std::string key_file = option_.key_directory + "/" + FileUtils::base_name(image_file) + ".key";
			if (!FileUtils::is_file(key_file)) {
				Logger::warn(title()) << "key file \'" << FileUtils::simple_name(key_file) << "\' doesn't exist" << std::endl;
				continue;
			}

			ImageData data;
			data.image_file = image_file;
			data.key_file = key_file;
			data.has_init_focal = true;
			data.init_focal = focal_length;
			data.image_ = nil;
			data.image_loaded = false;
			data.keys_loaded = false;
			data.camera.adjusted = false;
			image_data_.push_back(data);
		}

		// Create the match table
		match_table_ = MatchTable(num_of_images());
	}

	void SfM::load_match_table() {
		std::ifstream input(option_.match_table_file.c_str());
		if (input.fail()) {
			Logger::warn(title()) << "could not open match table file \'" << option_.match_table_file << "\'" << std::endl;
			return;
		}

		remove_all_matches();

		while (!input.eof()) {
			int img1, img2, num;
			input >> img1 >> img2 >> num;
			if (input.fail())
				break;

			set_match(img1, img2);

			/* Read the matches */
			std::vector<KeypointMatch> matches;
			for (int i = 0; i < num; ++i) {
				int k1, k2;
				input >> k1 >> k2;
				if (input.fail())
					break;

				KeypointMatch m;
				m.key_idx1 = k1;
				m.key_idx2 = k2;
				matches.push_back(m);
			}

			MatchIndex idx(img1, img2);
			match_table_.match_list(idx) = matches;
		}
	}

	void SfM::load_keys() {
		int num_images = num_of_images();
		for (int i = 0; i < num_images; i++) {
			image_data_[i].load_keys();
		}
	}


	/* Dump an output file containing information about the current state of the scene */
	void SfM::dump_output_file(const std::string& filename, 
		int num_images, int num_cameras, int num_points, 
		int* added_order, camera_params_t* cameras, vec3d* points, 
		vec3i* colors, std::vector<ImageKeyVector>& pt_views)
	{
		int num_visible_points = 0;
		for (int i = 0; i < num_points; i++) {
			if (pt_views[i].size() > 0)
				num_visible_points++;
		}
		if (num_visible_points == 0)
			return;

		std::ofstream output(filename.c_str());
		if (output.fail()) {
			Logger::warn(title()) << "could not create file \'" << filename << "\'" << std::endl;
			return;
		}

		/* Print version number */
		output << "# Bundle file v0.3" << std::endl;
		output << num_images << " " << num_visible_points << std::endl;

		/* Dump cameras */
		for (int i = 0; i < num_images; i++) {
			int idx = -1;
			for (int j = 0; j < num_cameras; j++) {
				if (added_order[j] == i) {
					idx = j;
					break;
				}
			}

			if (idx == -1) {
				output << "0 0 0" << std::endl;
				output << "0 0 0\n0 0 0\n0 0 0\n0 0 0\n";
			}
			else {
				output << cameras[idx].f << " " << cameras[idx].k[0] << " " << cameras[idx].k[1] << std::endl
					<< cameras[idx].R[0] << " " << cameras[idx].R[1] << " " << cameras[idx].R[2] << std::endl
					<< cameras[idx].R[3] << " " << cameras[idx].R[4] << " " << cameras[idx].R[5] << std::endl
					<< cameras[idx].R[6] << " " << cameras[idx].R[7] << " " << cameras[idx].R[8] << std::endl;

				double t[3];
				matrix_product(3, 3, 3, 1, cameras[idx].R, cameras[idx].t, t);
				matrix_scale(3, 1, t, -1.0, t);
				output << t[0] << " " << t[1] << " " << t[2] << std::endl;
			}
		}

		/* Dump points */
		for (int i = 0; i < num_points; i++) {
			int num_visible = (int)pt_views[i].size();

			if (num_visible > 0) {
				/* Position */
				output << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
				/* Color */
				output << colors[i].x << " " << colors[i].y << " " << colors[i].z << std::endl;

				int num_visible = (int)pt_views[i].size();
				output << num_visible << " ";
				for (int j = 0; j < num_visible; j++) {
					int img = added_order[pt_views[i][j].first];
					int key = pt_views[i][j].second;

					double x = image_data_[img].keys[key].x;
					double y = image_data_[img].keys[key].y;
					output << img << " " << key << " " << x << " " << y << " ";
				}
				output << std::endl;
			}
		}

#if 0
		/* Finally, dump all outliers */
		ImageKeyVector outliers;
		for (int i = 0; i < num_images; i++) {
			/* Find the index of this camera in the ordering */
			int idx = -1;
			for (int j = 0; j < num_cameras; j++) {
				if (added_order[j] == i) {
					idx = j;
					break;
				}
			}

			if (idx == -1) continue;

			int num_keys = GetNumKeys(i);
			for (int j = 0; j < num_keys; j++) {
				if (GetKey(i, j).m_extra == -2) {
					outliers.push_back(ImageKey(i, j));
				}
			}
		}

		int num_outliers = (int)outliers.size();
		fprintf(f, "%d\n", num_outliers);

		for (int i = 0; i < num_outliers; i++) {
			fprintf(f, "%d %d\n", outliers[i].first, outliers[i].second);
		}
#endif

	}



	/* Write point files to a ply file */
	void SfM::dump_points_to_ply(
		const std::string& filename, 
		int num_points, int num_cameras, 
		vec3d* points, vec3i* colors, camera_params_t* cameras)
	{
		int num_good_pts = 0;
		for (int i = 0; i < num_points; i++) {
			if (colors[i].x == 0x0 && colors[i].y == 0x0 && colors[i].z == 0xff)
				continue;
			++ num_good_pts;
		}
		if (num_good_pts == 0)
			return;

		std::ofstream output(filename.c_str());
		if (output.fail()) {
			Logger::warn(title()) << "could not create file \'" << filename << "\'" << std::endl;
			return;
		}

		Logger::out(title()) << "dumping points to ply file (" 
			<< num_good_pts << " points)" << std::endl;

		/* write the ply header */
		output << "ply" << std::endl
			<< "format ascii 1.0" << std::endl
			<< "element vertex " << num_good_pts + 2 * num_cameras << std::endl
			<< "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl
			<< "property uchar diffuse_red" << std::endl
			<< "property uchar diffuse_green" << std::endl
			<< "property uchar diffuse_blue" << std::endl
			<< "end_header" << std::endl;

		/* Now triangulate all the correspondences */
		for (int i = 0; i < num_points; i++) {
			if (colors[i].x == 0x0 && colors[i].y == 0x0 && colors[i].z == 0xff)
				continue;

			/* Output the vertex and its color*/
			output << points[i].x << " " 
				<< points[i].y << " " 
				<< points[i].z << " "
				<< iround(colors[i].x) << " "
				<< iround(colors[i].y) << " "
				<< iround(colors[i].z) << std::endl;
		}

		// export also the cameras
		for (int i = 0; i < num_cameras; i++) {
			double Rinv[9];
			matrix_invert(3, cameras[i].R, Rinv);

			double c[3];
			memcpy(c, cameras[i].t, 3 * sizeof(double));

			if ((i % 2) == 0)
				output << c[0] << " " << c[1] << " " << c[2] << " 0 255 0" << std::endl;
			else
				output << c[0] << " " << c[1] << " " << c[2] << " 255 0 0" << std::endl;

			double p_cam[3] = { 0.0, 0.0, -0.05 };
			double p[3];
			matrix_product(3, 3, 3, 1, Rinv, p_cam, p);

			p[0] += c[0];
			p[1] += c[1];
			p[2] += c[2];
			output << p[0] << " " << p[1] << " " << p[2] << " 255 255 0" << std::endl;
		}
	}


	void SfM::read_point_constraints(const std::string& file) {
		std::ifstream input(file.c_str());
		if (input.fail()) {
			Logger::warn(title()) << "could not create file \'" << file << "\'" << std::endl;
			return;
		}

		if (point_constraints_)
			delete point_constraints_;
		int num_points = (int)point_data_.size();
		point_constraints_ = new vec3d[num_points];

		std::string line;
		while (!input.eof()) {
			getline(input, line);
			if (input.fail())
				continue;

			if (line[0] == ' ' || line[0] == '%' || line[0] == '#')
				continue;  /* comment or whitespace */

			vec3d p0, p;
			std::istringstream string_stream(line);
			string_stream >> p0 >> p0;
			if (!string_stream.fail()) {
				int pt_idx = -1;
				double min_dist = DBL_MAX;
				for (int i = 0; i < num_points; i++) {
					double dsq = (point_data_[i].pos - p0).length2();
					if (dsq < min_dist) {
						pt_idx = i;
						min_dist = dsq;
					}
				}

				point_constraints_[pt_idx] = vec3d(p.x, p.y, -p.z);
			}
		}
	}


	void SfM::read_camera_constraints(const std::string& file) {
		std::ifstream input(file.c_str());
		if (input.fail()) {
			Logger::warn(title()) << "could not create file \'" << file << "\'" << std::endl;
			return;
		}

		std::string line;
		while (!input.eof()) {
			getline(input, line);
			if (input.fail())
				continue;

			if (line[0] == ' ' || line[0] == '%' || line[0] == '#')
				continue;  /* comment or whitespace */

			std::istringstream string_stream(line);
			int cam_idx;
			double x, y, z, xw, yw, zw;
			string_stream >> cam_idx >> x >> y >> z >> xw >> yw >> zw;

			if (x != -999.0) {
				image_data_[cam_idx].camera.constraints[0] = true;
				image_data_[cam_idx].camera.constraints[0] = x;
				image_data_[cam_idx].camera.constraints[0] = xw;
			}

			if (y != -999.0) {
				image_data_[cam_idx].camera.constraints[1] = true;
				image_data_[cam_idx].camera.constraints[1] = y;
				image_data_[cam_idx].camera.constraints[1] = yw;
			}

			if (z != -999.0) {
				image_data_[cam_idx].camera.constraints[2] = true;
				image_data_[cam_idx].camera.constraints[2] = z;
				image_data_[cam_idx].camera.constraints[2] = zw;
			}
		}
	}

}