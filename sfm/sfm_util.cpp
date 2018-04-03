#include "sfm_util.h"
#include <fstream>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <thread>

#include "../basic/basic_types.h"
#include "../basic/logger.h"
#include "../basic/file_utils.h"
#include "../image/image.h"
#include "../image/image_io.h"
#include "../math/matrix_driver.h"


namespace sfm {

	void read_bundle_file(const std::string& bundle_file,
		std::vector<camera_params_t> &cameras,
		std::vector<point_t> &points)
	{
		std::ifstream input(bundle_file.c_str());
		if (input.fail()) {
			std::cerr << "could not open file: \'" << bundle_file << "\'" << std::endl;
			return;
		}

		int num_images, num_points;

		std::string first_line;
		getline(input, first_line);
		std::istringstream line_stream(first_line);
		double bundle_version;
		if (first_line[0] == '#') {
			std::string dummy;
			char v;
			line_stream >> dummy >> dummy >> dummy >> v >> bundle_version;
			input >> num_images >> num_points;
		}
		else if (first_line[0] == 'v') {
			char v;
			line_stream >> v >> bundle_version;
			input >> num_images >> num_points;
		}
		else {
			bundle_version = 0.1;
			input >> num_images >> num_points;
		}

		std::cout << "[ReadBundleFile] Reading " << num_images << " images and " << num_points << " points" << std::endl;

		/* Read cameras */
		for (int i = 0; i < num_images; i++) {
			double focal_length, k0, k1;
			double R[9];
			double t[3];

			/* Focal lengths*/
			input >> focal_length >> k0 >> k1;
			/* Rotation */
			input >> R[0] >> R[1] >> R[2] >> R[3] >> R[4] >> R[5] >> R[6] >> R[7] >> R[8];
			/* Translation */
			input >> t[0] >> t[1] >> t[2];

			camera_params_t cam;
			cam.f = focal_length;
			cam.k[0] = k0;
			cam.k[1] = k1;
			memcpy(cam.R, R, sizeof(double) * 9);
			memcpy(cam.t, t, sizeof(double) * 3);

			/* Flip the scene if needed */
			if (bundle_version < 0.3) {
				R[2] = -R[2];
				R[5] = -R[5];
				R[6] = -R[6];
				R[7] = -R[7];
				t[2] = -t[2];
			}

			cameras.push_back(cam);
		}

		/* Read points */
		int total_num_visible = 0;
		for (int i = 0; i < num_points; i++) {
			point_t pt;

			/* Position */
			input >> pt.pos[0] >> pt.pos[1] >> pt.pos[2];
			/* Color */
			input >> pt.color[0] >> pt.color[1] >> pt.color[2];

			int num_visible;
			input >> num_visible;
			total_num_visible += num_visible;

			for (int j = 0; j < num_visible; j++) {
				int image, key;
				input >> image >> key;

				double x = 0, y = 0;
				if (bundle_version >= 0.3)
					input >> x >> y;

				view_t view;
				view.image = image;
				view.key = key;
				view.x = x;
				view.y = y;

				pt.views.push_back(view);
			}

			if (bundle_version < 0.3)
				pt.pos[2] = -pt.pos[2];

			if (num_visible > 0) {
				points.push_back(pt);
			}
		}
		
		std::cout << "Num visible: " << total_num_visible << std::endl;
	}

	void write_bundle_file(const std::string& bundle_file,
		const std::vector<camera_params_t> &cameras,
		const std::vector<point_t> &points)
	{
		std::ofstream output(bundle_file.c_str());
		if (output.fail()) {
			std::cerr << "could not open file: \'" << bundle_file << "\'" << std::endl;
			return;
		}

		int num_images = (int)cameras.size();
		int num_points = (int)points.size();

		/* Count the number of good images */
		int num_good_images = 0;
		int *map = new int[num_images];
		for (int i = 0; i < num_images; i++) {
			if (cameras[i].f == 0) {
				map[i] = -1;
				continue;
			}

			map[i] = num_good_images;
			num_good_images++;
		}

		std::cout << "[WriteBundleFile] Writing " << num_good_images << " images and " << num_points << " points..." << std::endl;

		output << "# Bundle file v0.3" << std::endl;
		output << num_good_images << " " << num_points << std::endl;

		/* Write cameras */
		for (int i = 0; i < num_images; i++) {
			if (cameras[i].f == 0)
				continue;

			/* Focal length */
			output << cameras[i].f << " 0.0 0.0" << std::endl;
			/* Rotation */
			const double *R = cameras[i].R;
			output << R[0] << " " << R[1] << " " << R[2] << " " 
				<< R[3] << " " << R[4] << " " << R[5] << " " 
				<< R[6] << " " << R[7] << " " << R[8] << std::endl;
			/* Translation */
			const double *t = cameras[i].t;
			output << t[0] << " " << t[1] << " " << t[2] << std::endl;
		}

		/* Write points */
		for (int i = 0; i < num_points; i++) {
			/* Position */
			const double *pos = points[i].pos;
			output << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;

			/* Color */
			const double *color = points[i].color;
			output << iround(color[0]) << " " << iround(color[1]) << " " << iround(color[2]) << std::endl;

			int num_visible = (int)points[i].views.size();
			output << num_visible << " ";

			for (int j = 0; j < num_visible; j++) {
				int view = map[points[i].views[j].image];
				assert(view >= 0 && view < num_good_images);
				int key = points[i].views[j].key;
				double x = points[i].views[j].x;
				double y = points[i].views[j].y;

				output << view << " " << key << " " << x << " " << y << " ";
			}

			output << std::endl;
		}
	}


	void read_list_file(const std::string& list_file, std::vector<std::string> &files)
	{
		std::ifstream input(list_file.c_str());
		if (input.fail()) {
			std::cerr << "could not open file: \'" << list_file << "\'" << std::endl;
			return;
		}

		std::string line;
		while (!input.eof()) {
			getline(input, line);
			std::istringstream line_stream(line);
			std::string file;
			line_stream >> file;
			if (!line_stream.fail())
				files.push_back(file);
		}
	}


	void write_pmvs(const std::string& pmvs_path,
		const std::string& txt_path,
		const std::string& list_file,
		const std::string& bundle_file,
		const std::string& pmvs_option_file,
		std::vector<std::string> images,
		std::vector<camera_params_t> &cameras)
	{
		int num_cameras = (int)cameras.size();
		int count = 0;
		for (int i = 0; i < num_cameras; i++) {
			if (cameras[i].f == 0.0)
				continue;

			std::ostringstream string_stream;
			string_stream << txt_path << "/" << std::setfill('0') << std::setw(4) << count << ".txt";
			std::string txt_file = string_stream.str();

			std::ofstream output(txt_file.c_str());
			if (output.fail()) {
				std::cerr << "could not open file: \'" << txt_file << "\'" << std::endl;
				return;
			}

			/* Compute the projection matrix */
			double focal = cameras[i].f;
			double *R = cameras[i].R;
			double *t = cameras[i].t;

			int w, h;
			ImageIO::query_image_size(images[i], w, h);

			double K[9] =
			{ -focal, 0.0, 0.5 * w - 0.5,
			0.0, focal, 0.5 * h - 0.5,
			0.0, 0.0, 1.0 };

			double Ptmp[12] =
			{ R[0], R[1], R[2], t[0],
			R[3], R[4], R[5], t[1],
			R[6], R[7], R[8], t[2] };

			double P[12];
			matrix_product(3, 3, 3, 4, K, Ptmp, P);
			matrix_scale(3, 4, P, -1.0, P);

			output << "CONTOUR" << std::endl;
			output << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << std::endl;
			output << P[4] << " " << P[5] << " " << P[6] << " " << P[7] << std::endl;
			output << P[8] << " " << P[9] << " " << P[10] << " " << P[11] << std::endl;
			count++;
		}

// 		std::cout << "@@ Sample command for running pmvs:\n"
// 			<< "   pmvs2 /PATH/pmvs_options.txt\n"
// 			<< "    - or - \n"
// 			<< "   use Dr. Yasutaka Furukawa's view clustering algorithm to generate a set of options files.\n"
// 			<< "       The clustering software is available at http://grail.cs.washington.edu/software/cmvs\n";

		// write the options file
		std::ofstream output(pmvs_option_file.c_str());
		if (!output.fail()) {
			output << "level 1" << std::endl
				<< "csize 2" << std::endl
				<< "threshold 0.7" << std::endl
				<< "wsize 7" << std::endl
				<< "minImageNum 3" << std::endl
				//<< "CPU 8" << std::endl // Liangliang: make full use of the CPU
				<< "CPU " << std::thread::hardware_concurrency() << std::endl
				<< "setEdge 0" << std::endl
				<< "useBound 0" << std::endl
				<< "useVisData 1" << std::endl
				<< "sequence -1" << std::endl
				<< "timages -1 0 " << count << std::endl
				<< "oimages -3" << std::endl;
		}
	}


	void write_vis_file(const std::string& vis_file,
		std::vector<camera_params_t> &cameras,
		std::vector<point_t> &points)
	{
		std::ofstream output(vis_file.c_str());
		if (output.fail()) {
			std::cerr << "could not open file: \'" << vis_file << "\'" << std::endl;
			return;
		}

		int nCameras = (int)cameras.size();
		int nPoints = (int)points.size();

		printf("Num cameras: %d\n", nCameras);

		/* Fill in the matches matrix */
		unsigned int *matches = new unsigned int[nCameras * nCameras];
		for (int i = 0; i < nCameras; i++) {
			for (int j = 0; j < nCameras; j++) {
				matches[i * nCameras + j] = 0;
			}
		}

		for (int i = 0; i < nPoints; i++) {
			// bool seen = false;
			int nViews = (int)points[i].views.size();
			for (int j = 0; j < nViews; j++) {
				int i1 = points[i].views[j].image;
				for (int k = j + 1; k < nViews; k++) {
					if (j == k) continue;
					int i2 = points[i].views[k].image;

					matches[i1 * nCameras + i2]++;
					matches[i2 * nCameras + i1]++;
				}
			}
		}

		output << "VISDATA" << std::endl;
		output << nCameras << std::endl;

		// write camera rows
		const unsigned int MATCH_THRESHOLD = 32;
		for (int i = 0; i < nCameras; i++) {
			std::vector<int> vis;
			for (int j = 0; j < nCameras; j++) {
				if (matches[i * nCameras + j] >= MATCH_THRESHOLD)
					vis.push_back(j);
			}

			int nVis = (int)vis.size();
			output << i << " " << nVis;

			for (int j = 0; j < nVis; j++) {
				output << " " << vis[j];
			}

			output << std::endl;
		}
	}


	void undistort_image(const std::string &in,
		const camera_params_t &camera,
		const std::string &out)
	{
		//printf("Undistorting image %s\n", in.c_str());
		Logger::out("SfM") << "Undistorting image " << FileUtils::simple_name(in) << std::endl;
		fflush(stdout);

		Image* img = ImageIO::read(in);
		if (!img) 
			return;
		//img_t *img = LoadJPEG(in.c_str());
		int w = img->width();
		int h = img->height();

		Image *img_out = new Image(img->color_encoding(), w, h);

		double f2_inv = 1.0 / (camera.f * camera.f);

		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				double x_c = x - 0.5 * w;
				double y_c = y - 0.5 * h;

				double r2 = (x_c * x_c + y_c * y_c) * f2_inv;
				double factor = 1.0 + camera.k[0] * r2 + camera.k[1] * r2 * r2;

				x_c *= factor;
				y_c *= factor;

				x_c += 0.5 * w;
				y_c += 0.5 * h;

				Colorf c;
				if (x_c >= 0 && x_c < w - 1 && y_c >= 0 && y_c < h - 1) {
					c = pixel_interpolate(img, x_c, y_c);
				}
				else {
					c = Colorf(0.0, 0.0, 0.0);
				}

				Memory::pointer ptr = img_out->pixel_base(x, y);
 				*(ptr) = iround(c.r());
 				*(ptr + 1) = iround(c.g());
 				*(ptr + 2) = iround(c.b());
			}
		}

		ImageIO::save(out, img_out);
		delete img;
		delete img_out;
	}

}
