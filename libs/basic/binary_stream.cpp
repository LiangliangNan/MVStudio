/*
* b_stream.cc
*/

#include "binary_stream.h"
#include "logger.h"
#include "assertions.h"

#include <cassert>
#include <iostream>
#include <cstdlib>



namespace IO {

	//____________________________________________________________________________

	void BinaryStream::detect_machine_endian() {

		Numeric::int32 data    = 0x04030201; 
		Numeric::uint8* pointer = (Numeric::uint8*)&data;

		if(
			pointer[0] == 1 && 
			pointer[1] == 2 &&
			pointer[2] == 3 &&
			pointer[3] == 4
			) {
				machine_endian_ = IO_LITTLE_ENDIAN;
				//            Logger::out("BinaryStream") 
				//                << "machine type: little endian" << std::endl;
		} else if(
			pointer[0] == 4 && 
			pointer[1] == 3 &&
			pointer[2] == 2 &&
			pointer[3] == 1) {
				machine_endian_ = IO_BIG_ENDIAN;
				//            Logger::out("BinaryStream") 
				//                << "machine type: big endian" << std::endl;
		} else {
			Logger::err("BinaryStream") 
				<< "invalid processor type" << std::endl;
			// Should not occur, nowadays machines are either little or
			// big endian.
			::abort() ;
		}
	}

	void BinaryStream::init_convert_tables(void) {
		if(stream_endian_ != machine_endian_) {
			endian_convert_64_[0] = 7;
			endian_convert_64_[1] = 6;
			endian_convert_64_[2] = 5;
			endian_convert_64_[3] = 4;	     
			endian_convert_64_[4] = 3;
			endian_convert_64_[5] = 2;
			endian_convert_64_[6] = 1;
			endian_convert_64_[7] = 0;

			endian_convert_32_[0] = 3;
			endian_convert_32_[1] = 2;
			endian_convert_32_[2] = 1;
			endian_convert_32_[3] = 0;

			endian_convert_16_[0] = 1;
			endian_convert_16_[1] = 0;
		} else {
			endian_convert_64_[0] = 0;
			endian_convert_64_[1] = 1;
			endian_convert_64_[2] = 2;
			endian_convert_64_[3] = 3;	     
			endian_convert_64_[4] = 4;
			endian_convert_64_[5] = 5;
			endian_convert_64_[6] = 6;
			endian_convert_64_[7] = 7;	     

			endian_convert_32_[0] = 0;
			endian_convert_32_[1] = 1;
			endian_convert_32_[2] = 2;
			endian_convert_32_[3] = 3;

			endian_convert_16_[0] = 0;
			endian_convert_16_[1] = 1;	     
		}
	}


	//____________________________________________________________________________

	BinaryInputStream::BinaryInputStream(
		const std::string& file_name, int stream_endian_in
		) : BinaryStream(stream_endian_in) {
			record_OK_    = true ;
			record_count_ = 0;
			input_ = new std::ifstream(
				file_name.c_str(),
				std::fstream::in | std::fstream::binary
				);
			owns_input_ = true ;
			if(!input_) {
				delete input_; input_ = nil;
			}
			has_record_markers_ = true ;
	}

	BinaryInputStream::BinaryInputStream(
		std::istream& input, int stream_endian_in
		) : BinaryStream(stream_endian_in) {
			record_OK_    = true ;
			record_count_ = 0;
			input_ = &input ;
			owns_input_ = false ;
			has_record_markers_ = true ;
	}

	bool BinaryInputStream::OK(void) const {
		return (record_OK_ && (input_ != nil)) ; 
	}

	bool BinaryInputStream::more(void) const {
		return !input_->eof();
	}

	void BinaryInputStream::begin_record(void) {
		if(input_->eof()) {
			record_OK_ = false ;
		} else {
			if(has_record_markers_) {
				(*this) >> count1_;
			}
		}
	}

	void BinaryInputStream::end_record(void){
		record_count_++;
		if(input_->eof()) {
			record_OK_ = false ;
		} else {
			if(has_record_markers_) {
				(*this) >> count2_;
				if(count1_ != count2_) {
					record_OK_ = false ;
					Logger::err("BinaryStream")
						<< "Invalid record length in record #"
						<< record_count_ << std::endl;
					Logger::err("BinaryStream") 
						<< " count1=" << count1_ 
						<< " count2=" << count2_ << std::endl;
				}
			}
		}
	}

	BinaryInputStream::~BinaryInputStream(void) {
		if(input_ != nil && owns_input_) {
			delete input_ ;
		}
		input_ = nil ;
	}

	BinaryInputStream& 
		BinaryInputStream::read_array(Numeric::uint8* data, int nbr) {
			for(int i=0; i<nbr; i++)
				(*this) >> data[i];
			return *this;
	}

	BinaryInputStream& 
		BinaryInputStream::read_array(Numeric::int16* data, int nbr) {
			for(int i=0; i<nbr; i++)
				(*this) >> data[i];
			return *this;	
	}

	BinaryInputStream& 
		BinaryInputStream::read_array(Numeric::int32* data, int nbr) {
			for(int i=0; i<nbr; i++)
				(*this) >> data[i];
			return *this;	
	}

	BinaryInputStream& 
		BinaryInputStream::read_array(Numeric::float32* data, int nbr) {
			for(int i=0; i<nbr; i++)
				(*this) >> data[i];
			return *this;	
	}

	BinaryInputStream& 
		BinaryInputStream::read_array(Numeric::float64* data, int nbr) {
			for(int i=0; i<nbr; i++)
				(*this) >> data[i];
			return *this;	
	}

	BinaryInputStream&
		BinaryInputStream::read_record(Numeric::uint8* data, int nbr) {
			begin_record();
			read_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryInputStream&
		BinaryInputStream::read_record(Numeric::int16* data, int nbr) {
			begin_record();
			read_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryInputStream&
		BinaryInputStream::read_record(Numeric::int32* data, int nbr) {
			begin_record();
			read_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryInputStream&
		BinaryInputStream::read_record(Numeric::float32* data, int nbr) {
			begin_record();
			read_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryInputStream&
		BinaryInputStream::read_record(Numeric::float64* data, int nbr) {
			begin_record();
			read_array(data, nbr);
			end_record();
			return *this;
	}
	//____________________________________________________________________________

	BinaryOutputStream::BinaryOutputStream(
		const std::string& file_name, int stream_endian_in
		) : BinaryStream(stream_endian_in) {
			output_ = new std::ofstream(
				file_name.c_str(),
				std::fstream::out | std::fstream::trunc | std::fstream::binary
				) ;
			if(!output_) {
				delete output_; output_ = nil ;
			}
			owns_output_ = true ;
	}

	BinaryOutputStream::BinaryOutputStream(
		std::ostream& output, int stream_endian_in
		) : BinaryStream(stream_endian_in) {
			output_ = &output ;
			owns_output_ = false ;
	}


	bool BinaryOutputStream::OK(void) const {
		return (output_ != nil);
	}

	void BinaryOutputStream::begin_record(void) {
		if(has_record_markers_) {
			count_ = 0;
			pos_ = output_->tellp() ;
			write_marker(count_);
		}
	}

	void BinaryOutputStream::end_record(void)  {
		if(has_record_markers_) {
			write_marker(count_);            
			Numeric::int32 pos = output_->tellp();   
			output_->seekp(pos_) ;
			write_marker(count_);            
			output_->seekp(pos) ;
		}
	}

	BinaryOutputStream::~BinaryOutputStream(void) {
		if(output_ != nil && owns_output_) {
			delete output_ ;
		}
		output_ = nil ;
	}

	void BinaryOutputStream::write_marker(Numeric::int32 x) {
		ogf_assert(output_ != nil);
		char* x_orig = (char *)&x;
		for(unsigned int i=0; i<sizeof(x); i++) {
			output_->put(x_orig[endian_convert_32_[i]])  ;
		}
	}

	BinaryOutputStream& BinaryOutputStream::write_array(
		Numeric::uint8* data, int nbr
		) {
			for(int i=0; i<nbr; i++)
				(*this) << data[i];
			return *this;
	}

	BinaryOutputStream& BinaryOutputStream::write_array(
		Numeric::int16* data, int nbr
		) {
			for(int i=0; i<nbr; i++)
				(*this) << data[i];
			return *this;	
	}

	BinaryOutputStream& BinaryOutputStream::write_array(
		Numeric::int32* data, int nbr
		) {
			for(int i=0; i<nbr; i++)
				(*this) << data[i];
			return *this;	
	}

	BinaryOutputStream& BinaryOutputStream::write_array(
		Numeric::float32* data, int nbr
		) {
			for(int i=0; i<nbr; i++)
				(*this) << data[i];
			return *this;	
	}

	BinaryOutputStream& BinaryOutputStream::write_array(
		Numeric::float64* data, int nbr
		) {
			for(int i=0; i<nbr; i++)
				(*this) << data[i];
			return *this;	
	}

	BinaryOutputStream& BinaryOutputStream::write_record(
		Numeric::uint8* data, int nbr
		) {
			begin_record();
			write_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryOutputStream& BinaryOutputStream::write_record(
		Numeric::int16* data, int nbr
		) {
			begin_record();
			write_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryOutputStream& BinaryOutputStream::write_record(
		Numeric::int32* data, int nbr
		) {
			begin_record();
			write_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryOutputStream& BinaryOutputStream::write_record(
		Numeric::float32* data, int nbr
		) {
			begin_record();
			write_array(data, nbr);
			end_record();
			return *this;
	}

	BinaryOutputStream& BinaryOutputStream::write_record(
		Numeric::float64* data, int nbr
		) {
			begin_record();
			write_array(data, nbr);
			end_record();
			return *this;
	}

	//__________________________________________________________
}

