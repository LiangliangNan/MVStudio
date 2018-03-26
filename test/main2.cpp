#include <iostream>
#include "../basic/logger.h"
#include "../basic/file_utils.h"
#include "../basic/basic_types.h"
#include "../3rd_numeric/SuiteSparse-4.2.1/SuiteSparse_config/SuiteSparse_config.h"


int main(int argc, char* argv[])
{
    size_t a, b;
    b = 1;
    b = a + 1;

    std::cout << "sizeof(size_t): " << sizeof(size_t) * 8<< std::endl;

    std::cout << "sizeof(int): " << sizeof(int) * 8 << std::endl;
	std::cout << "sizeof(long): " << sizeof(long) * 8 << std::endl;
	std::cout << "sizeof(long int): " << sizeof(long int) * 8 << std::endl;
	std::cout << "sizeof(SuiteSparse_long): " << sizeof(SuiteSparse_long) * 8 << std::endl;
	
	std::cout << "sizeof(int8): " << sizeof(Numeric::int8) * 8 << std::endl;
	std::cout << "sizeof(int16): " << sizeof(Numeric::int16) * 8 << std::endl;
	std::cout << "sizeof(int32): " << sizeof(Numeric::int32) * 8 << std::endl;
	std::cout << "sizeof(int64): " << sizeof(Numeric::int64) * 8 << std::endl;
	std::cout << "sizeof(uint8): " << sizeof(Numeric::uint8) * 8 << std::endl;
	std::cout << "sizeof(uint16): " << sizeof(Numeric::uint16) * 8 << std::endl;
	std::cout << "sizeof(uint32): " << sizeof(Numeric::uint32) * 8 << std::endl;
	std::cout << "sizeof(uint64): " << sizeof(Numeric::uint64) * 8 << std::endl;
	std::cout << "sizeof(float32): " << sizeof(Numeric::float32) * 8 << std::endl;
	std::cout << "sizeof(float64): " << sizeof(Numeric::float64) * 8 << std::endl;







    Logger::initialize();
    Logger::instance()->set_value(Logger::LOG_REGISTER_FEATURES, "*"); // log everything
    Logger::instance()->set_value(Logger::LOG_FILE_NAME, "MeshStudio.log");


    //FileUtils::set_current_working_directory("d:");

    //std::cout << "current working directory: " << FileUtils::get_current_working_directory() << std::endl;

    //std::string name = "./test_dir";
    //if (FileUtils::is_directory(name))
    //	std::cout << "\'" << name << "\' is a directory" << std::endl;
    //if (FileUtils::is_file(name))
    //	std::cout << "\'" << name << "\' is a file" << std::endl;

    //std::vector<std::string> files_recursive, files_no_recursive;
    //FileUtils::get_files(name, files_recursive, true) ;
    //FileUtils::get_files(name, files_no_recursive, false) ;

    //std::vector<std::string> sub_dirs_recursive, sub_dirs_no_recursive;
    //FileUtils::get_subdirectories(name, sub_dirs_recursive, true) ;
    //FileUtils::get_subdirectories(name, sub_dirs_no_recursive, false) ;



    //name = "./bug.obj";
    //if (FileUtils::is_directory(name))
    //	std::cout << "\'" << name << "\' is a directory" << std::endl;
    //if (FileUtils::is_file(name))
    //	std::cout << "\'" << name << "\' is a file" << std::endl;


    return 0;
}

