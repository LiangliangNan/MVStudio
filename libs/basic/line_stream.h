
#ifndef _LINE_STREAM_H_
#define _LINE_STREAM_H_

#include "basic_types.h"
#include "assertions.h"

#include <iostream>
#include <sstream>



namespace IO {

	//____________________________________________________________________________

	class LineInputStream {
	public:
		LineInputStream(std::istream& in) : in_(in), line_in_(nil) {   }
		~LineInputStream() {
			delete line_in_ ; 
			line_in_ = nil ;
		}

		bool eof() const { return in_.eof() ; }
		bool eol() const { return line_in_ == nil || line_in_->eof() ; }
		//bool ok() const { return in_ != nil; } // Liangliang: changed for vs2013
		bool ok() const { return !in_.fail(); }

		void get_line() {
			getline(in_, buffer_);
			//in_.getline(buffer_, 65536) ;
			delete line_in_ ; 
			line_in_ = new std::istringstream(buffer_) ;
		}

		std::istream& line() { 
			ogf_assert(line_in_ != nil) ;
			return *line_in_ ; 
		}

		const std::string& current_line() const {
			return buffer_;
		}

		template <class T> LineInputStream& operator>>(T& param) {
			*line_in_ >> param;
			return *this;
		}

	private:
		std::istream& in_ ;
		std::istringstream* line_in_ ;
		std::string buffer_;
		//char buffer_[65536] ;
	} ;
}


typedef IO::LineInputStream     LineInputStream;



#endif
