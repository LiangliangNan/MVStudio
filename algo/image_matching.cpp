
#include "image_matching.h"

#include "project.h"

#include "../basic/file_utils.h"
#include "../basic/progress.h"
#include "../basic/logger.h"
#include "../basic/stop_watch.h"

//include the header files for SiftGPU and SiftMatchGPU
#include "../3rd_party/SiftGPU/SiftGPU.h"
#include "../3rd_party/SiftGPU/SiftMatch.h"

#include <algorithm>




ImageMatching::ImageMatching(Project* proj)
	: project_(proj)
{
}

ImageMatching::~ImageMatching() {
}


struct ImageFeature {
	std::vector<SiftGPU::SiftKeypoint>	keys;
	std::vector<float>			descriptors;
};

std::vector<ImageFeature>  image_features;



void ImageMatching::apply() {
	if (!project_ || !project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return;
	}

	image_features.clear();
	StopWatch w;
	Logger::out(title()) << "extracting key points..." << std::endl;
	extract_key_points();
	Logger::out(title()) << "done. time: " << w.elapsed() << std::endl;

	w.start();
	Logger::out(title()) << "matching the images..." << std::endl;
	match_key_points();
	Logger::out(title()) << "done. time: " << w.elapsed() << std::endl;
}



void ImageMatching::extract_key_points() {
	if (!project_ || !project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return;
	}

	std::ofstream output(project_->sfm_list_file.c_str());
	if (output.fail()) {
		Logger::err(title()) << "could not write focal length file \'" << project_->sfm_list_file << "\'" << std::endl;
		return;
	}

	SiftGPU* sift = CreateNewSiftGPU();

	char * argv[] = { "-fo", "-1", "-v", "0", "-tc2", "7680", "-b", "-nomc" };//
	//-fo -1    staring from -1 octave 
	//-v 1      only print out # feature and overall time
	//-tc <num> set a soft limit to number of detected features
	// NEW: by default SiftGPU will try to fit the cap of GPU memory, and reduce the working 
	// dimension so as to not allocate too much. This feature can be disabled by -nomc
	int argc = sizeof(argv)/sizeof(char*);
	sift->ParseParam(argc, argv);
	//sift->SetMaxDimension(8000);

	///////////////////////////////////////////////////////////////////////
	//Only the following parameters can be changed after initialization (by calling ParseParam). 
	//-dw, -ofix, -ofix-not, -fo, -unn, -maxd, -b
	//to change other parameters at runtime, you need to first unload the dynamically loaded library
	//reload the library, then create a new siftgpu instance

	image_features.clear();

	ProgressLogger progress(int(project_->images.size()));
	for (std::size_t i = 0; i < project_->images.size(); ++i) {
		if (progress.is_canceled()) {
			image_features.clear();
			delete sift;
			return;
		}
		if (project_->images[i].ignored)
			continue;

		//Create a context for computation, and SiftGPU will be initialized automatically 
		//The same context can be used by SiftMatchGPU.
		// When you mix your own OpenGL code with SiftGPU, you need to recall 
		// CreateContextGL before calling SiftGPU functions, which will implicitly activate
		// the internal OpenGL context.
		if (sift->CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED) {
			Logger::err(title()) << "SiftGPU is not supported}" << std::endl;
			break;
		}

		const std::string& name = project_->images[i].file;
		if (!sift->RunSIFT(name.c_str())) {
			Logger::warn(title()) << "processing image \'" << FileUtils::simple_name(name) << "\' failed" << std::endl;
			project_->set_ignore_image(name, true);
			continue;
		}

		std::string sift_file = project_->sfm_keys_dir + '/' + FileUtils::base_name(name) + ".key";
		sift->SaveSIFT(sift_file.c_str());

		int image_width, image_height;
		sift->GetImageDimension(image_width, image_height);
		float focal_length_in_pixels = 1.2f * std::max(image_width, image_height);
		output << FileUtils::get_absolute_path(project_->images[i].file) 
			<< " w: " << image_width 
			<< " h: " << image_height << " f: " << focal_length_in_pixels << std::endl;

		int num = sift->GetFeatureNum();	// get feature count
		Logger::out(title()) << FileUtils::simple_name(name) << ": " << num << " key points" << std::endl;

		std::vector<SiftGPU::SiftKeypoint> keys;
		std::vector<float> descriptors;

		if (num > 0) {
			keys.resize(num);
			descriptors.resize(num * 128);
			//reading back feature vectors is faster than writing files
			//if you don't need keys or descriptors, just put NULLs here
			sift->GetFeatureVector(&keys[0], &descriptors[0]);
			//this can be used to write your own sift file.  
		}

		ImageFeature feature;
		feature.keys = keys;
		feature.descriptors = descriptors;
		image_features.push_back(feature);

		progress.next();
	}

	delete sift;
}


void ImageMatching::match_key_points() {
	if (!project_ || !project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return;
	}

	// save the matches into a file
	std::string match_table_file = project_->sfm_match_table_file;
	std::ofstream output(match_table_file.c_str());
	if (output.fail()) {
		Logger::err(title()) << "could not write matches file \'" << match_table_file << "\'" << std::endl;
		return;
	}

	SiftMatchGPU* matcher = new SiftMatchGPU;

	//Before initialization, you can choose between GLSL, and CUDA(if compiled). 
	//matcher->SetLanguage(SiftMatchGPU::SIFTMATCH_CUDA); // +i for the (i+1)-th device

	//Set descriptors to match, the first argument must be either 0 or 1
	//if you want to use more than 4096 or less than 4096
	//call matcher->SetMaxSift() to change the limit before calling setdescriptor

	int num = (int)image_features.size();
	ProgressLogger progress(num * (num - 1) / 2);

	int count = 0;
	for ( int i = 0; i < num; ++i) {
		const std::vector<SiftGPU::SiftKeypoint>& keys1 = image_features[i].keys;
		const std::vector<float>& descriptors1 = image_features[i].descriptors;
		int num1 = (int)keys1.size();
		if (num1 == 0)
			continue;
		
		for ( int j = i+1; j < num; ++j) {
			const std::vector<SiftGPU::SiftKeypoint>& keys2 = image_features[j].keys;
			const std::vector<float>& descriptors2 = image_features[j].descriptors;
			int num2 = (int)keys2.size();
			if (num2 == 0)
				continue;

			if (progress.is_canceled()) {
				delete matcher;
				return;
			}

			// If you already have an OpenGL Context, call matcher->VerifyContextGL() instead
			matcher->CreateContextGL();
			//Set descriptors to match, the first argument must be either 0 or 1
			//if you want to use more than 4096 or less than 4096
			//call matcher->SetMaxSift() to change the limit before calling setdescriptor
			matcher->SetDescriptors(0, num1, &descriptors1[0]); //image 1
			matcher->SetDescriptors(1, num2, &descriptors2[0]); //image 2

			//match and get result.    
			int(*match_buf)[2] = new int[num1][2];
			//use the default thresholds. Check the declaration in SiftGPU.h
			int num_match = matcher->GetSiftMatch(num1, match_buf);
			output << i << "   " << j << std::endl << num_match << std::endl;

			//enumerate all the feature matches
			for (int k = 0; k < num_match; ++k) {
				//How to get the feature matches: 

				// get the index in the sift keys of the image
				int id1 = match_buf[k][0];
				int id2 = match_buf[k][1];
				output << id1 << " " << id2 << " " << std::endl;

				//key1 in the first image matches with key2 in the second image
// 				const SiftGPU::SiftKeypoint & key1 = keys1[id1];
// 				const SiftGPU::SiftKeypoint & key2 = keys2[id2];
// 				(key1.x, key1.y);
// 				(key2.x, key2.y);
					
			}
			output << std::endl;

			delete[] match_buf;

			progress.next();
		}
	}

	image_features.clear();
	delete matcher;
}
