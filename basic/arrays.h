
#ifndef _ARRAYS_H_
#define _ARRAYS_H_


#include "basic_types.h"
#include "assertions.h"

#define OGF_ARRAY_BOUND_CHECK

#ifdef OGF_ARRAY_BOUND_CHECK
#define ogf_array_check(index, size) ogf_range_assert(index, 0, size-1) ;
#else 
#define ogf_array_check(index, size) 
#endif


template <class T> class Array1d {
public:
	typedef Array1d<T> thisclass ;

	Array1d(int size = 0, int alignment = 1) {
		data_ = nil ;
		base_mem_ = nil ;
		size_ = 0 ;
		alignment_ = 1 ;
		allocate(size, alignment) ;
	} 

	inline ~Array1d() { deallocate() ; }

	/** does not preserve previous values stored in the array */
	void allocate(int size, int alignment = 1) {
		deallocate() ;
		if(size != 0) {
			base_mem_ = (Memory::pointer)malloc(size * sizeof(T) + alignment -1) ;
			Memory::pointer p = base_mem_ ;
			while(Numeric::uint64(p) % alignment) {  p++ ; }
			data_ = (T*)p ;
			for(int i=0; i<size; i++) {
				// Direct call to the constructor, see dlist.h for more explanations.
				new(&data_[i])T() ;                    
			}
		}
		size_ = size ;
		alignment_ = alignment ;
	}

	void set_all(const T& value) {
		for(unsigned int i=0; i<size_; i++) {
			data_[i] = value ;
		}
	}

	T& operator()(unsigned int index) {
		ogf_array_check(index, size_) ;
		return data_[index] ;
	}

	const T& operator()(unsigned int index) const {
		ogf_array_check(index, size_) ;
		return data_[index] ;
	}

	T& operator[](unsigned int index) {
		ogf_array_check(index, size_) ;
		return data_[index] ;
	}

	const T& operator[](unsigned int index) const {
		ogf_array_check(index, size_) ;
		return data_[index] ;
	}

	T& from_linear_index(unsigned int index) {
		ogf_array_check(index, size_) ;
		return data_[index] ;
	}

	const T& from_linear_index(unsigned int index) const {
		ogf_array_check(index, size_) ;
		return data_[index] ;
	}

	unsigned int size() const { return size_ ; }
	unsigned int alignment() const { return alignment_ ; }

	void clear() { allocate(0) ; }

	/** low-level access, for experts only. */
	const T* data() const { return data_ ; }

	/** low-level access, for experts only. */
	T* data() { return data_ ; }

	unsigned int mem_usage() const {
		return size_ * sizeof(T) + sizeof(thisclass) ;
	}

protected:
	T* data_ ;
	unsigned int size_ ;
	Memory::pointer base_mem_ ;
	unsigned int alignment_ ;

protected:
	void deallocate() {
		if(size_ != 0) {
			for(unsigned int i=0; i<size_; i++) {
				// direct call to the destructor, see dlist.h for more explanations
				data_[i].~T() ;
			}
			free(base_mem_) ;
			data_ = nil ;
			base_mem_ = nil ;
			size_ = 0 ;
			alignment_ = 1 ;
		}
	}

private:
	Array1d(const thisclass& rhs) ;
	thisclass& operator=(const thisclass& rhs) ;
} ;

#ifdef OGF_OS_WINDOWS
// Fix for MSVC
template class BASIC_API Array1d<double>;
#endif

//_________________________________________________________

template <class T> class Array2d {
public:
	typedef Array2d<T> thisclass ;

	Array2d(int size1 = 0, int size2 = 0) {
		data_ = nil ;
		allocate(size1, size2) ;
	} 

	inline ~Array2d() { delete[] data_ ; data_ = nil ; }

	/** does not preserve previous values stored in the array */
	void allocate(int size1, int size2) {
		delete[] data_ ; 
		int size = size1 * size2 ;
		data_ = (size == 0) ? nil : new T[size]; 
		size_[0] = size1 ;
		size_[1] = size2 ;
	}

	void set_all(const T& value) {
		int size = size_[0] * size_[1] ;
		for(int i=0; i<size; i++) {
			data_[i] = value ;
		}
	}

	T& operator()(int index0, int index1) {
		ogf_array_check(index0, size_[0]) ;
		ogf_array_check(index1, size_[1]) ;
		return data_[index1 * size_[0] + index0] ;
	}

	const T& operator()(int index0, int index1) const {
		ogf_array_check(index0, size_[0]) ;
		ogf_array_check(index1, size_[1]) ;
		return data_[index1 * size_[0] + index0] ;
	}

	T& from_linear_index(int index) {
		ogf_array_check(index, size_[0] * size_[1]) ;
		return data_[index] ;
	}

	const T& from_linear_index(int index) const {
		ogf_array_check(index, size_[0] * size_[1]) ;
		return data_[index] ;
	}

	int size(int dim) const { 
		ogf_array_check(dim, 2) ;
		return size_[dim] ;
	}

	void clear() { allocate(0,0) ; }

	const T* data() const { return data_ ; }
	T* data() { return data_ ; }

private:
	T* data_ ;
	int size_[2] ;

private:
	Array2d(const thisclass& rhs) ;
	thisclass& operator=(const thisclass& rhs) ;
} ;


//_________________________________________________________


template <class T> class Array3d {
public:
	typedef Array3d<T> thisclass ;

	Array3d(int size1 = 0, int size2 = 0, int size3 = 0) {
		data_ = nil ;
		allocate(size1, size2, size3) ;
	} 

	inline ~Array3d() { delete[] data_ ; data_ = nil ; }

	/** does not preserve previous values stored in the array */
	void allocate(int size1, int size2, int size3) {
		delete[] data_ ; 
		int size = size1 * size2 * size3 ;
		data_ = (size == 0) ? nil : new T[size]; 
		size_[0] = size1 ;
		size_[1] = size2 ;
		size_[2] = size3 ;
	}

	void set_all(const T& value) {
		int size = size_[0] * size_[1] * size_[2] ;
		for(int i=0; i<size; i++) {
			data_[i] = value ;
		}
	}

	T& operator()(int index0, int index1, int index2) {
		ogf_array_check(index0, size_[0]) ;
		ogf_array_check(index1, size_[1]) ;
		ogf_array_check(index2, size_[2]) ;
		return data_[index0 + size_[0] * (index1 + size_[1] * index2)] ;
	}

	const T& operator()(int index0, int index1, int index2) const {
		ogf_array_check(index0, size_[0]) ;
		ogf_array_check(index1, size_[1]) ;
		ogf_array_check(index2, size_[2]) ;
		return data_[index0 + size_[0] * (index1 + size_[1] * index2)] ;
	}

	T& from_linear_index(int index) {
		ogf_array_check(index, size_[0] * size_[1] * size_[2]) ;
		return data_[index] ;
	}

	const T& from_linear_index(int index) const {
		ogf_array_check(index, size_[0] * size_[1] * size_[2]) ;
		return data_[index] ;
	}

	int size(int dim) const { 
		ogf_array_check(dim, 3) ;
		return size_[dim] ;
	}

	void clear() { allocate(0,0,0) ; }

	const T* data() const { return data_ ; }
	T* data() { return data_ ; }

private:
	T* data_ ;
	int size_[3] ;

private:
	Array3d(const thisclass& rhs) ;
	thisclass& operator=(const thisclass& rhs) ;
} ;

//_________________________________________________________

/**
* Fixed size array with bound checking, copy and copy constructor.
* FixedArray1d can also be used as an Attribute, in contrast with
* regular C++ arrays.
*/
template <class T, int N> class FixedArray1d {
public:
	enum { dimension = N } ;
	typedef FixedArray1d<T,N> thisclass ;

	FixedArray1d() { }
	FixedArray1d(const thisclass& rhs) { copy(rhs) ; }
	thisclass& operator=(const thisclass& rhs) { copy(rhs) ; return *this ; }

	const T& operator()(int i) const { 
		ogf_array_check(i, N) ;
		return data_[i] ; 
	}

	T& operator()(int i) { 
		ogf_array_check(i, N) ;
		return data_[i] ; 
	}        

	void set_all(const T& value) {
		for(int i=0; i<N; i++) {
			data_[i] = value ;
		}
	}

	const T* data() const { return data_ ; }
	T* data() { return data_ ; }

protected:
	void copy(const thisclass& rhs) {
		for(int i=0; i<N; i++) {
			data_[i] = rhs.data_[i] ;
		}
	}

private:
	T data_[N] ;
} ;



#endif

