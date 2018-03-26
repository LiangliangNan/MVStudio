
#ifndef __BASIC_OS_TEXT_UTILS__
#define __BASIC_OS_TEXT_UTILS__

#include "basic_common.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>



namespace TextUtils {

	class BASIC_API Environment {
	public:
		bool has_variable(const std::string& name) const ;
		std::string value(const std::string& name) const ;
		const std::vector<std::string>& values(const std::string& name) const ;
		void set_value(const std::string& name, const std::string& value) ;
		void set_values(const std::string& name, const std::vector<std::string>& values) ;
		void append_value(const std::string& name, const std::string& value) ;
		void append_values(const std::string& name, const std::vector<std::string>& values) ;
		void clear_value(const std::string& name) ;
		void clear() ;
		void print(std::ostream& out) const ;

	private:
		typedef std::map< std::string, std::vector<std::string> > EnvMap ;        
		EnvMap data_ ;
	} ;

	inline std::ostream& operator<<(std::ostream& out, const Environment& env) {	
		env.print(out) ;
		return out ;
	}

	void BASIC_API read_environment_file(
		const std::string& file_name,
		Environment& environment
		) ;

	void BASIC_API find_and_replace(
		std::istream& in, std::ostream& out,
		const Environment& env
		) ;

	void BASIC_API concatenate(
		std::istream& in, std::ostream& out
		) ;


	void BASIC_API flip_slashes(std::string& s) ;

}


#endif
