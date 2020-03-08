
#ifndef _BINARY_STREAM_H_
#define _BINARY_STREAM_H_

#include "basic_common.h"
#include "basic_types.h"

#include <iostream>
#include <fstream>
#include <string>


namespace IO {

	//____________________________________________________________________________

    const int IO_LITTLE_ENDIAN = 0 ;
    const int IO_BIG_ENDIAN    = 1 ;

	/**
	* A characteristic of processors named endian corresponds to the
	* way they store numbers. A little endian processor stores
	* less significant bytes first, and a big endian processor
	* stores most significant bytes first.
	* Most Unix machines are big endian, except those having
	* a DEC Alpha processor. Intel processors are little endian.
	* This causes portability problems for binary files to be
	* exchanged between machines of the two types. The classes
	* BinaryInputStream and BinaryOutputStream make it possible
	* to deal with this problem. \\
	* We have not used XDR because XDR assumes that the external
	* format is always big endian, which would have prevented us
	* to read files from PC softwares (such as 3DS binary files). 
	*/

	class BASIC_API BinaryStream  {

	public:
		BinaryStream(int stream_endian_in = IO_BIG_ENDIAN) ;

		void stream_endian(int stream_endian_in) ;
		int stream_endian() const ;
		int machine_endian() const ;

		/**
		* checks the flag positioned by set_had_record_markers().
		*/
		bool has_record_markers() const ;

		/** 
		* some FORTRAN files have their records surrounded by two
		* integers and some not. This function enables BinaryStream
		* to read such FORTRAN files. Default value is true.  
		*/
		void set_has_record_markers(bool b) ;

	protected:
		void detect_machine_endian() ;
		void init_convert_tables() ;

		/** endian of the processor */
		int     machine_endian_;
		/** endian of the file */
		int     stream_endian_;
		/** conversion table for 64 bits words */
		Numeric::uint8 endian_convert_64_[8];
		/** conversion table for 32 bits words */    
		Numeric::uint8 endian_convert_32_[4];
		/** conversion table for 16 bits words */    
		Numeric::uint8 endian_convert_16_[2]; 

		bool has_record_markers_ ;
	} ;

	//____________________________________________________________________________

	/**
	* This class enables binary files to be read, while taking
	* into account endian problems (see BinaryStream).
	*/

	class BASIC_API BinaryInputStream : public BinaryStream {

	public:
		BinaryInputStream(
			const std::string& file_name, int endian_in = IO_BIG_ENDIAN
			);

		BinaryInputStream(
			std::istream& input, int endian_in = IO_BIG_ENDIAN
			);

		~BinaryInputStream();

		BinaryInputStream& operator>>(Numeric::uint8&);
		BinaryInputStream& operator>>(Numeric::int16&);
		BinaryInputStream& operator>>(Numeric::int32&);
		BinaryInputStream& operator>>(Numeric::uint16&);
		BinaryInputStream& operator>>(Numeric::uint32&);
		BinaryInputStream& operator>>(Numeric::float32&);
		BinaryInputStream& operator>>(Numeric::float64&); 

		/** Reads a data bloc without converting */    
		BinaryInputStream& read_opaque_data(void*, int size); 
		BinaryInputStream& read_opaque_data(void*, int size, int nbr);

		/** checks the state of the stream */
		bool OK(void) const; 

		/** returns true if data remains to be read */    
		bool more(void) const; 

		/** 
		* FORTRAN data files are structured into records,
		* bounded by two integers indicating the size of
		* the record. These two functions enable these integers
		* to be read, and to use them as a validity check. If
		* they differ, subsequent calls to OK() return false.
		* Note that if set_had_record_markers() has been called
		* with false, the records are supposed to be continuously
		* written in the file (without markers).
		*/

		void begin_record(void);
		void end_record(void);

		/** array read functions */

		BinaryInputStream& read_array(
			Numeric::uint8*   data, int nbr
			); 
		BinaryInputStream& read_array(
			Numeric::int16*  data, int nbr
			); 
		BinaryInputStream& read_array(
			Numeric::int32*  data, int nbr
			);
		BinaryInputStream& read_array(
			Numeric::float32*  data, int nbr
			);
		BinaryInputStream& read_array(
			Numeric::float64* data, int nbr
			);

		/** reading a record made of an array */

		BinaryInputStream& read_record(
			Numeric::uint8*   data, int nbr
			); 
		BinaryInputStream& read_record(
			Numeric::int16*  data, int nbr
			); 
		BinaryInputStream& read_record(
			Numeric::int32*  data, int nbr
			); 
		BinaryInputStream& read_record(
			Numeric::float32*  data, int nbr
			);
		BinaryInputStream& read_record(
			Numeric::float64* data, int nbr
			);

		unsigned long tell() const {
			return input_->tellg() ;
		}

		void seek(unsigned long p) {
			input_->seekg(p) ;
		}

		void seek(unsigned long p, std::ios_base::seekdir dir) {
			input_->seekg(p, dir) ;
		}


	private:
		std::istream* input_ ;
		bool owns_input_ ;
		/** No pb encountered when reading the last record */
		bool record_OK_;
		/** Sentry integer preceding the current record */
		Numeric::int32 count1_;
		/** Sentry integer following the current record */
		Numeric::int32 count2_;
		/** Number of read records */
		Numeric::int32 record_count_; 
	};

	//____________________________________________________________________________

	/**
	* Enables binary files to be written, while taking
	* into account endian problems. 
	*/

	class BASIC_API BinaryOutputStream : public BinaryStream {

	public:
		BinaryOutputStream(
			const std::string& file_name, int stream_endian = IO_BIG_ENDIAN
			);
		BinaryOutputStream(
			std::ostream& output, int stream_endian = IO_BIG_ENDIAN
			);

		~BinaryOutputStream();

		BinaryOutputStream& operator<<(Numeric::uint8);
		BinaryOutputStream& operator<<(Numeric::int16);
		BinaryOutputStream& operator<<(Numeric::int32);
		BinaryOutputStream& operator<<(Numeric::uint16);
		BinaryOutputStream& operator<<(Numeric::uint32);
		BinaryOutputStream& operator<<(Numeric::float32);
		BinaryOutputStream& operator<<(Numeric::float64);	

		/** Write a data bloc without conversion */
		BinaryOutputStream& write_opaque_data(void*, int size);
		BinaryOutputStream& write_opaque_data(void*, int size, int nbr);

		/** checks the state of the stream */
		bool OK(void) const;   

		/** 
		* FORTRAN data files are structured into records,
		* bounded by two integers indicating the size of
		* the record. These two functions enable these integers
		* to be read, and to use them as a validity check. If
		* they differ, subsequent calls to OK() return false.
		* Note that if set_had_record_markers() has been called
		* with false, the records are supposed to be continuously
		* written in the file (without markers).
		*/
		void begin_record(void);          
		void end_record(void);

		/** Writing arrays */
		BinaryOutputStream& write_array(
			Numeric::uint8*   data, int nbr
			);
		BinaryOutputStream& write_array(
			Numeric::int16*  data, int nbr
			); 
		BinaryOutputStream& write_array(
			Numeric::int32*  data, int nbr
			);
		BinaryOutputStream& write_array(
			Numeric::float32*  data, int nbr
			);
		BinaryOutputStream& write_array(
			Numeric::float64* data, int nbr
			);	

		/** Writing a record made of an array */
		BinaryOutputStream& write_record(
			Numeric::uint8*   data, int nbr
			); 
		BinaryOutputStream& write_record(
			Numeric::int16*  data, int nbr
			); 
		BinaryOutputStream& write_record(
			Numeric::int32*  data, int nbr
			); 
		BinaryOutputStream& write_record(
			Numeric::float32*  data, int nbr
			);
		BinaryOutputStream& write_record(
			Numeric::float64* data, int nbr
			);

	private:
		/** Write a sentry integer (@see InputStream) */
		void    write_marker(Numeric::int32 value);
		std::ostream*   output_ ;
		bool owns_output_ ;
		/** Size of the current record */
		Numeric::int32 count_;
		/**
		* Position of the current record relative to the
		* beginning of the file. This is used to write count_
		* at that location when end_record is called.
		*/
		Numeric::int32 pos_; 
	};

	//____________________________________________________________________________

	inline BinaryStream::BinaryStream(int stream_endian_in) {
		stream_endian_ = stream_endian_in ;
		has_record_markers_ = true ;
		detect_machine_endian() ;
		init_convert_tables() ;
	}

	inline int BinaryStream::machine_endian() const {
		return machine_endian_ ;
	}

	inline int BinaryStream::stream_endian() const {
		return stream_endian_;
	}

	inline void BinaryStream::stream_endian(int endian_in) {
		stream_endian_ = endian_in ;
		init_convert_tables() ;
	}

	inline bool BinaryStream::has_record_markers() const {
		return has_record_markers_ ;
	}

	inline void BinaryStream::set_has_record_markers(bool b) {
		has_record_markers_ = b ;
	}

	//____________________________________________________________________________

	inline BinaryInputStream& BinaryInputStream::operator>>(
		Numeric::uint8& x
		) {
			char* dest = (char *)&x; 
			//        input_->get(*dest) ;
			input_->read(dest, 1) ;
			return *this ;
	}

	inline BinaryInputStream& BinaryInputStream::operator>>(
		Numeric::int16& x
		) {
			char* dest = (char *)&x; 
			for(unsigned int i=0; i<sizeof(x); i++) {
				input_->read(&dest[endian_convert_16_[i]], 1) ;
			}
			return *this ;
	}

	inline BinaryInputStream& BinaryInputStream::operator>>(
		Numeric::int32& x
		) {
			char* dest = (char *)&x; 
			for(unsigned int i=0; i<sizeof(x); i++) {
				input_->read(&dest[endian_convert_32_[i]], 1) ;
			}
			return *this ;
	}

	inline BinaryInputStream& BinaryInputStream::operator>>(
		Numeric::uint16& x
		) {
			char* dest = (char *)&x; 
			for(unsigned int i=0; i<sizeof(x); i++) {
				input_->read(&dest[endian_convert_16_[i]], 1) ;
			}
			return *this ;
	}

	inline BinaryInputStream& BinaryInputStream::operator>>(
		Numeric::uint32& x
		) {
			char* dest = (char *)&x; 
			for(unsigned int i=0; i<sizeof(x); i++) {
				input_->read(&dest[endian_convert_32_[i]], 1) ;
			}
			return *this ;
	}


	inline BinaryInputStream& BinaryInputStream::operator>>(
		Numeric::float32& x
		) {
			char* dest = (char *)&x; 
			for(unsigned int i=0; i<sizeof(x); i++) {
				input_->read(&dest[endian_convert_32_[i]], 1) ;
			}
			return *this ;
	}

	inline BinaryInputStream& BinaryInputStream::operator>>(
		Numeric::float64& x
		) {
			char* dest = (char *)&x; 
			for(unsigned int i=0; i<sizeof(x); i++) {
				input_->read(&dest[endian_convert_64_[i]], 1) ;
			}
			return *this ;
	}

	inline BinaryInputStream& BinaryInputStream::read_opaque_data(
		void* p, int size
		) {
			char* dest = (char *)p ;
			input_->read(dest, size) ;
			return *this ;
	} 

	inline BinaryInputStream& BinaryInputStream::read_opaque_data(
		void* ptr, int size, int nbr
		) {
			return read_opaque_data(ptr, size * nbr);
	}

	//____________________________________________________________________________


	inline BinaryOutputStream& BinaryOutputStream::operator<<(
		Numeric::uint8 x
		) {
			char* src = (char*)&x ;
			output_->put(*src) ;
			count_ += sizeof(x) ;
			return *this ;
	}

	inline BinaryOutputStream& BinaryOutputStream::operator<<(
		Numeric::int16 x
		) {
			char* src = (char*)&x ;
			for(unsigned int i=0; i<sizeof(x); i++) {
				output_->write(&src[endian_convert_16_[i]],1) ;
			}
			count_ += sizeof(x) ;
			return *this ;
	}

	inline BinaryOutputStream& BinaryOutputStream::operator<<(
		Numeric::int32 x
		) {
			char* src = (char*)&x ;
			for(unsigned int i=0; i<sizeof(x); i++) {
				output_->write(&src[endian_convert_32_[i]],1) ;
			}
			count_ += sizeof(x) ;
			return *this ;
	}

	inline BinaryOutputStream& BinaryOutputStream::operator<<(
		Numeric::uint16 x
		) {
			char* src = (char*)&x ;
			for(unsigned int i=0; i<sizeof(x); i++) {
				output_->write(&src[endian_convert_16_[i]],1) ;
			}
			count_ += sizeof(x) ;
			return *this ;
	}

	inline BinaryOutputStream& BinaryOutputStream::operator<<(
		Numeric::uint32 x
		) {
			char* src = (char*)&x ;
			for(unsigned int i=0; i<sizeof(x); i++) {
				output_->write(&src[endian_convert_32_[i]],1) ;
			}
			count_ += sizeof(x) ;
			return *this ;
	}

	inline BinaryOutputStream& BinaryOutputStream::operator<<(
		Numeric::float32 x
		) {
			char* src = (char*)&x ;
			for(unsigned int i=0; i<sizeof(x); i++) {
				output_->write(&src[endian_convert_32_[i]],1) ;
			}
			count_ += sizeof(x) ;
			return *this ;
	}

	inline BinaryOutputStream& BinaryOutputStream::operator<<(
		Numeric::float64 x
		) {
			char* src = (char*)&x ;
			for(unsigned int i=0; i<sizeof(x); i++) {
				output_->write(&src[endian_convert_64_[i]],1) ;
			}
			count_ += sizeof(x) ;
			return *this ;
	}	

	inline BinaryOutputStream& BinaryOutputStream::write_opaque_data(
		void* p, int size
		) {
			char* src = (char*)p ;
			output_->write(src, size) ;
			count_ += size ;
			return *this ;
	}

	inline BinaryOutputStream& 
		BinaryOutputStream::write_opaque_data(void* ptr, int size, int nbr) {
			return write_opaque_data(ptr, size * nbr);
	}

	//____________________________________________________________________________
}



#endif
