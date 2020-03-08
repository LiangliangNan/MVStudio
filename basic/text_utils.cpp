
#include "text_utils.h"
#include "assertions.h"
#include "basic_types.h"

#include <fstream>
#include <iostream>

#ifdef WIN32
#include <string.h>
#else
#include <strings.h>
#endif


namespace TextUtils {

	// this function avoids the \n or c10 or c13 errors at end of line
	static void getline(std::istream& in, char* out, int maxsize) {
		::memset(out,0,maxsize) ;
		char c_temp = 0 ;
		int i = 0 ;
		in.get(c_temp) ;
		while((c_temp == 10 || c_temp == 13) && !in.eof()) {
			in.get(c_temp) ;
		}

		while(c_temp != 10 && c_temp != 13 && !in.eof()) {
			out[i] = c_temp ;
			i++ ;
			ogf_assert(i < maxsize) ;
			in.get(c_temp) ;
		}
		out[i] = '\0' ;
	}


	bool Environment::has_variable(const std::string& name) const {
		EnvMap::const_iterator it = data_.find(name) ;
		return it != data_.end() ;
	}

	std::string Environment::value(const std::string& name) const {
		const std::vector<std::string>& vals = values(name) ;
		std::string result ;
		for(unsigned int i=0; i<vals.size(); i++) {
			if(i != 0) {
				result += " " ;
			}
			result += vals[i] ;
		}
		return result ;
	}

	const std::vector<std::string>& Environment::values(const std::string& name) const {
		EnvMap::const_iterator it = data_.find(name) ;
		ogf_assert(it != data_.end()) ;
		return it->second ;
	}

	void Environment::set_value(const std::string& name, const std::string& value) {
		data_[name].clear() ;
		append_value(name,value) ;
	}

	void Environment::set_values(
		const std::string& name, const std::vector<std::string>& values
		) {
			data_[name].clear() ;
			append_values(name,values) ;
	}

	void Environment::append_value(const std::string& name, const std::string& value) {
		if(value[0] == '$') {
			std::string var_name = value.substr(1,value.length() - 1) ;
			// Fabien 21/11/2003
			if (has_variable(var_name)) {
				// it's a local var
				append_values(name,values(var_name)) ;
			}
			else {
				// it's a windows var like (QTDIR)
				data_[name].push_back(value) ;
			}
			//ogf_assert(has_variable(var_name)) ;
			//append_values(name,values(var_name)) ;
		} else {
			data_[name].push_back(value) ;
		}
	}

	void Environment::append_values(
		const std::string& name, const std::vector<std::string>& values
		) {
			for(unsigned int i=0; i<values.size(); i++) {
				append_value(name,values[i]) ;
			}
	}

	void Environment::clear_value(const std::string& name) {
		ogf_assert(has_variable(name)) ;
		data_[name].clear() ;
	}

	void Environment::clear() {
		data_.clear() ;
	}

	void Environment::print(std::ostream& out) const {
		for(
			Environment::EnvMap::const_iterator it = data_.begin();
			it != data_.end(); it++
			) {
				// Fabien Boutantin: 2004/03/08 
				// added double quotes to output
				out << it->first << "= \" " ;
				for(unsigned int i=0; i<it->second.size(); i++) {
					out << it->second[i] << " " ;
				}
				out << "\"" << std::endl ;
		}
	}


	//_________________________________________________________________________________________

	static void read_environment_variable(
		const std::string& variable,
		Environment& environment
		) {
			char buff[1024] ;
            strcpy(buff, variable.c_str()) ;

			char* p_equal = strchr(buff, '=') ;
			ogf_assert(p_equal != nil) ;
			*p_equal = '\0' ;

			std::string variable_name(buff) ;

			char* p_value = p_equal + 1 ;
			p_value = strchr(p_value, '"') ;
			ogf_assert(p_value != nil) ;
			p_value++ ;

			char* p_value_end = strchr(p_value, '"') ;
			ogf_assert(p_value_end != nil) ;
			*p_value_end = '\0' ;

			std::string variable_values_str(p_value) ;
			std::vector<std::string> variable_values ;
			::String::split_string(variable_values_str, ' ', variable_values) ;
			environment.set_values(variable_name, variable_values) ;
	}

	void read_environment_file(
		const std::string& file_name,
		Environment& environment
		) {
			std::ifstream in(file_name.c_str()) ;
			ogf_assert(in) ;

			char buff[1024] ;
			std::string line ;
			bool append = false ;

			while(in) {
				getline(in,buff,1024) ;

				char* p = strchr(buff, '#') ;
				if(p != nil) {
					*p = '\0' ;
				}

				if(buff[0] == '\0') {
					continue ;
				}

				if(buff[strlen(buff) - 1] == '\\') {
					append = true ;
					buff[strlen(buff) - 1] = '\0' ;
				} else {
					append = false ;
				}

				line += std::string(buff) ;

				if(!append) {
					read_environment_variable(line, environment) ;
					line.clear() ;
				}

			}
	}


	void find_and_replace(
		std::istream& in, std::ostream& out,
		const Environment& env
		) {
			char buff[1024] ;
			std::string var_name ;
			while(in) {
				getline(in,buff,1024) ;
				std::string in_string(buff) ;
				std::string out_string ;
				bool in_var = false ;
				for(unsigned int i=0; i<in_string.length(); i++) {
					char cur_char = in_string[i] ;
					if(in_var) {
						if(cur_char == '%') {
							in_var = false ;
							if(env.has_variable(var_name)) {
								out_string += env.value(var_name) ;
							} else {
								// Maybe a variable for MSDOS, pass through...
								out_string += ("%" + var_name + "%") ;
							}
							var_name.clear() ;
						} else {
							var_name += ::String::char_to_string(cur_char) ;
						}
					} else {
						if(cur_char == '%') {
							in_var = true ;
						} else {
							out_string += ::String::char_to_string(cur_char) ;
						}
					}
				}
				out << out_string << std::endl ;
			}
	}

	void concatenate(
		std::istream& in, std::ostream& out
		) {
			char buff[1024] ;
			while(in) {
				getline(in,buff,1024) ;
				out << buff << std::endl ;
			}
	}

	void flip_slashes(std::string& s) {
		for(unsigned int i=0; i<s.length(); i++) {
			if(s[i] == '\\') {
				s[i] = '/' ;
			}
		}
	}

}
