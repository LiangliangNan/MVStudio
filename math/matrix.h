

#ifndef __MATH_LINEAR_ALGEBRA_MATRIX__
#define __MATH_LINEAR_ALGEBRA_MATRIX__

#include "../basic/assertions.h"
#include <iostream>



/**
* A class for representing matrices where the coefficients
* may be of arbitrary types.
*/

template <class FT, int DIM> class Matrix ;
template <class FT, int DIM> Matrix<FT, DIM> operator*(
	FT op1, const Matrix<FT, DIM>& op2
	) ;


template <class FT, int DIM> 
class Matrix 
{
public:
	Matrix() ;
	void load_zero() ;
	void load_identity() ;

	FT& operator()(int i, int j) ;
	const FT& operator()(int i, int j) const ;

	Matrix<FT, DIM>& operator+=(const Matrix<FT, DIM>& rhs) ;
	Matrix<FT, DIM>& operator-=(const Matrix<FT, DIM>& rhs) ;
	Matrix<FT, DIM>& operator*=(FT rhs) ;
	Matrix<FT, DIM>& operator/=(FT rhs) ;

	Matrix<FT, DIM> operator+(const Matrix<FT, DIM>& op2) const ;
	Matrix<FT, DIM> operator-(const Matrix<FT, DIM>& op2) const ;
	Matrix<FT, DIM> operator*(const Matrix<FT, DIM>& op2) const ;

	Matrix<FT, DIM> inverse() const ;
	Matrix<FT, DIM> transpose() const ;

	// cast to Scalar array
	operator FT*() { return &(coeff_[0][0]); }

	// cast to const Scalar array
	operator const FT*() const { return &(coeff_[0][0]); }

	// Routines for interfacing with Fortran, OpenGL etc...
	const FT* data() const { return &(coeff_[0][0]) ; }
	FT* data() { return &(coeff_[0][0]) ; }

	void get_lower_triangle(FT* store) {
		for(unsigned int i=0; i<DIM; i++) {
			for(unsigned int j=0; j<=i; j++) {
				*store++ = coeff_[i][j] ;
			}
		}
	}

private:
	FT coeff_[DIM][DIM] ;
} ;

//_______________________________________________________________________

template <class FT, int DIM> inline 
Matrix<FT, DIM>::Matrix() {
	load_identity() ;
}


template <class FT, int DIM> inline 
FT& Matrix<FT, DIM>::operator()(int i, int j) {
	ogf_assert(i >= 0 && i < DIM) ;
	ogf_assert(j >= 0 && j < DIM) ;
	return coeff_[i][j] ;
}

template <class FT, int DIM> inline 
const FT& Matrix<FT, DIM>::operator()(int i, int j) const {
	ogf_assert(i >= 0 && i < DIM) ;
	ogf_assert(j >= 0 && j < DIM) ;
	return coeff_[i][j] ;
}

template <class FT, int DIM> inline 
Matrix<FT, DIM>& Matrix<FT, DIM>::operator+=(const Matrix<FT, DIM>& rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] += rhs.coeff_[i][j] ;
		}
	}
	return *this ;
}

template <class FT, int DIM> inline 
Matrix<FT, DIM>& Matrix<FT, DIM>::operator-=(const Matrix<FT, DIM>& rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] -= rhs.coeff_[i][j] ;
		}
	}
	return *this ;
}

template <class FT, int DIM> inline 
Matrix<FT, DIM>& Matrix<FT, DIM>::operator*=(FT rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] *= rhs ;
		}
	}
	return *this ;
}

template <class FT, int DIM> inline 
Matrix<FT, DIM>& Matrix<FT, DIM>::operator/=(FT rhs) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] /= rhs ;
		}
	}
	return *this ;
}

template <class FT, int DIM> inline 
Matrix<FT, DIM> Matrix<FT, DIM>::operator+(const Matrix<FT, DIM>& op2) const {
	Matrix<FT, DIM> result = *this ;
	result += op2 ;
	return result ;
}

template <class FT, int DIM> inline 
Matrix<FT, DIM> Matrix<FT, DIM>::operator-(const Matrix<FT, DIM>& op2) const {
	Matrix<FT, DIM> result = *this ;
	result -= op2 ;
	return result ;

}

template <class FT, int DIM> inline 
Matrix<FT, DIM> Matrix<FT, DIM>::operator*(const Matrix<FT, DIM>& op2) const {
	Matrix<FT, DIM> result ;
	result.load_zero() ;
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			for(int k=0; k<DIM; k++) {
				result.coeff_[i][j] += coeff_[i][k] * op2.coeff_[k][j] ;
			}
		}
	}
	return result ;
}

template <class FT, int DIM> inline 
std::ostream& operator << (std::ostream& output, const Matrix<FT, DIM>& m) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			output << m(i,j) << " " ;
		}
	}    
	return output ;
}

template <class FT, int DIM> inline 
std::istream& operator >> (std::istream& input, Matrix<FT, DIM>& m) {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			input >> m(i,j) ;
		}
	}    
	return input ;
}



//_______________________________________________________________________

template <class FT, int DIM> void
Matrix<FT, DIM>::load_zero() {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] = FT(0) ;
		}
	}
}

template <class FT, int DIM> void
Matrix<FT, DIM>::load_identity() {
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			coeff_[i][j] = (i==j) ? FT(1) : FT(0) ;
		}
	}
}

template <class FT, int DIM> Matrix<FT, DIM> 
Matrix<FT, DIM>::inverse() const {
	FT val, val2;
	int i, j, k, ind;
	Matrix<FT, DIM> tmp = (*this);
	Matrix<FT, DIM> result;

	result.load_identity();

	for (i = 0; i != DIM; i++) {
		val = tmp(i,i);			/* find pivot */
		ind = i;
		for (j = i + 1; j != DIM; j++) {
			if (std::fabs(tmp(j,i)) > std::fabs(val)) {
				ind = j;
				val = tmp(j,i);
			}
		}

		if (ind != i) {			
			for (j = 0; j != DIM; j++) {
				val2 = result(i,j);
				result(i,j) = result(ind,j);
				result(ind,j) = val2;           /* swap columns */
				val2 = tmp(i,j);
				tmp(i,j) = tmp(ind,j);
				tmp(ind,j) = val2;
			}
		}

		ogf_assert(val != 0.0);

		for (j = 0; j != DIM; j++) {
			tmp(i,j)    /= val;
			result(i,j) /= val;
		}

		for (j = 0; j != DIM; j++) {		
			if (j == i)
				continue;                       /* eliminate column */
			val = tmp(j,i);
			for (k = 0; k != DIM; k++) {
				tmp(j,k)     -= tmp(i,k)     * val;
				result(j,k)  -= result(i,k)  * val;
			}
		}
	}

	return result;
}

template <class FT, int DIM> Matrix<FT, DIM> 
Matrix<FT, DIM>::transpose() const {
	Matrix<FT, DIM> result ;
	for(int i=0; i<DIM; i++) {
		for(int j=0; j<DIM; j++) {
			result(i,j) = (*this)(j,i) ;
		}
	}
	return result ;
}

template <class FT, int N> inline
void mult(const Matrix<FT,N>& M, const FT* x, FT* y) {
	for(int i=0; i<N; i++) {
		y[i] = 0 ;
		for(int j=0; j<N; j++) {
			y[i] += M(i,j)*x[j] ;
		}
	}
}


#endif
