
#ifndef _BASIC_SMART_POINTER_H_
#define _BASIC_SMART_POINTER_H_

#include "assertions.h"


/**
* Automatic memory management using reference counting. 
* This class can be used with classes inheriting
* the Counted class.
* @see Counted
*/

template <class T> class SmartPointer {

public:
	SmartPointer() ;
	SmartPointer(T* ptr) ;
	SmartPointer(const SmartPointer<T>& rhs) ;
	~SmartPointer() ;

	SmartPointer<T>& operator=(T* ptr) ;
	SmartPointer<T>& operator=(const SmartPointer<T>& rhs) ;

	T* get() { return pointer_; }
	const T* get() const { return pointer_; }

	/**
	* Makes the current instance a nil pointer ("forgets" the
	* current reference). 'p.forget();' is a shorthand for
	* 'p = nil;'
	*/
	void forget() ;

	T* operator->() const ;
	T& operator*() const ;
	operator T*() const ;

	/**
	* 'if(p.is_nil()) {...}' is a shorthand for 'if(p == nil) {...}'
	*/
	bool is_nil() const ;

protected:
private:
	T* pointer_ ;
} ;

//____________________________________________________________________________

template <class T> inline 
SmartPointer<T>::SmartPointer() : pointer_(nil) {
}

template <class T> inline 
SmartPointer<T>::SmartPointer(T* ptr) : pointer_(ptr) {
	T::ref(pointer_) ;
}

template <class T> inline 
SmartPointer<T>::SmartPointer(
							  const SmartPointer<T>& rhs
							  ) : pointer_(rhs) {
								  T::ref(pointer_) ;
} 

template <class T> inline
SmartPointer<T>::~SmartPointer() {
	T::unref(pointer_) ;
}

template <class T> inline
SmartPointer<T>& SmartPointer<T>::operator=(T* ptr) {
	if(ptr != pointer_) {
		T::unref(pointer_) ;
		pointer_ = ptr ;
		T::ref(pointer_) ;
	}
	return *this ;
}

template <class T> inline
SmartPointer<T>& SmartPointer<T>::operator=(const SmartPointer<T>& rhs) {
	T* rhs_p = rhs ; 
	if(rhs_p != pointer_) {
		T::unref(pointer_) ;
		pointer_ = rhs_p ;
		T::ref(pointer_) ;
	}
	return *this ;
}

template <class T> inline
void SmartPointer<T>::forget() {
	T::unref(pointer_) ;
	pointer_ = nil ;
}

template <class T> inline
T* SmartPointer<T>::operator->() const {
	ogf_assert(pointer_ != nil) ;
	return pointer_ ;
}

template <class T> inline
T& SmartPointer<T>::operator*() const {
	ogf_assert(pointer_ != nil) ;
	return *pointer_ ;
}

template <class T> inline 
SmartPointer<T>::operator T*() const {
	return pointer_ ;
}

template <class T> inline 
bool SmartPointer<T>::is_nil() const {
	return (pointer_ == nil) ;
}




#endif
