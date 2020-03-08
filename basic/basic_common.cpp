#include "basic_common.h"
#include "logger.h"
#include "progress.h"
#include "assertions.h"

namespace mpl {

	static class init_lib_basic
	{
	public:
		init_lib_basic()
		{
			Logger::instance()->set_value(Logger::LOG_REGISTER_FEATURES, "*"); // log everything
			 //Logger::instance()->set_value(Logger::LOG_FILE_NAME, log_file_);
			 //Logger::instance()->set_value("log_features",
			 //	"EigenSolver;MeshBuilder;MeshParameterizer\
			 //	LinearSolver");
		}

		~init_lib_basic() {
			Logger::terminate();
			Progress::terminate();
		}

	}  init;
}

