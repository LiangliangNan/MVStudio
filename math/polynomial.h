
#ifndef _MATH_SYMBOLIC_POLYNOMIAL__
#define _MATH_SYMBOLIC_POLYNOMIAL__


#include "../basic/assertions.h"
#include "factor.h"

#include <cmath>
#include <algorithm>


// N degree polynomial 

template< int N >
class Polynomial {
public:
	enum { degree = N };
	enum { nb_coeffs = N + 1 };

	Polynomial(void);
	Polynomial(double coeff[]);

	template<int N2>
	Polynomial(const Polynomial<N2>& P);

	void	set_coefficients(double coeff[]);
	double* coefficients();
	const double* coefficients() const;

	void	set_coefficient(int i, double c);
	double& coefficient(int i);
	const double& coefficient(int i) const;

	double operator()(double x) const;
	double eval(double x) const;

	double integral(double x_min, double x_max) const;

	bool operator == (const Polynomial& p) const;
	bool operator != (const Polynomial& p) const;

	bool is_zero(void) const;
	void set_zero(void);

	void print(void) const;

	template<int N2>
	Polynomial& operator  = (const Polynomial<N2> &p);
	Polynomial& operator += (const Polynomial& p);
	Polynomial& operator -= (const Polynomial& p);
	Polynomial  operator -  (void) const;
	Polynomial  operator +  (const Polynomial& p) const;
	Polynomial  operator -  (const Polynomial& p) const;
	template<int N2>
	Polynomial<N + N2>  operator *  (const Polynomial<N2>& p) const;

	Polynomial& operator += (double s);
	Polynomial& operator -= (double s);
	Polynomial& operator *= (double s);
	Polynomial& operator /= (double s);
	Polynomial  operator +  (double s) const;
	Polynomial  operator -  (double s) const;
	Polynomial  operator *  (double s) const;
	Polynomial  operator /  (double s) const;

	Polynomial scale(double s) const;
	Polynomial shift(double t) const;

	Polynomial<N - 1> derivative(void) const;
	Polynomial<N + 1> integral(void) const;

	Polynomial& add_scaled(const Polynomial& p, double scale);

	static void negate(const Polynomial& in, Polynomial& out);
	static void subtract(const Polynomial& p1, const Polynomial& p2, Polynomial& q);
	static void scale(const Polynomial& p, double w, Polynomial& q);
	static void add_scaled(const Polynomial& p1, double w1, const Polynomial& p2, double w2, Polynomial& q);
	static void add_scaled(const Polynomial& p1, const Polynomial& p2, double w2, Polynomial& q);
	static void add_scaled(const Polynomial& p1, double w1, const Polynomial& p2, Polynomial& q);

	void get_solutions(double value, std::vector<double>& roots, double eps) const;
	int  get_solutions(double value, double* roots, double eps) const; // returns the num of roots

	static Polynomial BSpline_component(int i);

private:
	double coeff_[nb_coeffs];
};




//////////////////////////////////////////////////////////////////////////

template<int N>
Polynomial<N>::Polynomial(void) {
	memset(coeff_, 0, sizeof(double) * nb_coeffs);
}


template<int N>
Polynomial<N>::Polynomial(double coeff[]) {
	for (int i = 0; i < nb_coeffs; ++i) {
		coeff_[i] = coeff[i];
	}
}


template<int N>
template<int N2>
Polynomial<N>::Polynomial(const Polynomial<N2>& P) {
	memset(coeff_, 0, sizeof(double) * nb_coeffs);
	for (int i = 0; i <= N && i <= N2; i++) {
		coeff_[i] = P.coefficient(i);
	}
}


template<int N>
void Polynomial<N>::set_coefficients(double coeff[]) {
	for (int i = 0; i < nb_coeffs; ++i) {
		coeff_[i] = coeff[i];
	}
}


template<int N>
double* Polynomial<N>::coefficients() {
	return coeff_;
}


template<int N>
const double* Polynomial<N>::coefficients() const {
	return coeff_;
}


template<int N>
void Polynomial<N>::set_coefficient(int i, double c) {
	assert(i >= 0 && i < nb_coeffs);
	coeff_[i] = c;
}


template<int N>
double& Polynomial<N>::coefficient(int i) {
	return coeff_[i];
}


template<int N>
const double& Polynomial<N>::coefficient(int i) const {
	return coeff_[i];
}


template<int N>
template<int N2>
Polynomial<N>& Polynomial<N>::operator  = (const Polynomial<N2> &p) {
	int d = N < N2 ? N : N2;
	memset(coeff_, 0, sizeof(double) * nb_coeffs);
	memcpy(coeff_, p.coefficients(), sizeof(double)*(d + 1));
	return *this;
}


template<int N>
Polynomial<N - 1> Polynomial<N>::derivative(void) const{
	Polynomial<N - 1> p;
	for (int i = 0; i < N; i++) {
		p.coefficient(i) = coeff_[i + 1] * (i + 1);
	}
	return p;
}


template<int N>
Polynomial<N + 1> Polynomial<N>::integral(void) const{
	Polynomial<N + 1> p;
	p.coefficients()[0] = 0;
	for (int i = 0; i <= N; i++) {
		p.coefficients()[i + 1] = coeff_[i] / (i + 1); 
	}
	return p;
}


template<> double Polynomial< 0 >::operator() (double t) const { return coeff_[0]; }

template<> double Polynomial< 1 >::operator() (double t) const { return coeff_[0] + coeff_[1] * t; }

template<> double Polynomial< 2 >::operator() (double t) const { return coeff_[0] + (coeff_[1] + coeff_[2] * t) * t; }


template<int N>
double Polynomial<N>::operator() (double x) const{
	double result = coeff_[N];
	for (int i = N - 1; i >= 0; i--) 
		result = result * x + coeff_[i];
	return result;
}

template<int N>
double Polynomial<N>::eval(double x) const {
	double result = coeff_[N];
	for (int i = N - 1; i >= 0; i--) {
		result = x * result + coeff_[i];
	}
	return result;
}

template<int N>
double Polynomial<N>::integral(double x_min, double x_max) const
{
	double result = 0;
	double x1 = x_min;
	double x2 = x_max;
	for (int i = 0; i <= N; i++){
		result += coeff_[i] * (x2 - x1) / (i + 1);
		if (x1 > -DBL_MAX && x1 < DBL_MAX)
			x1 *= x_min;
		if (x2 > -DBL_MAX && x2 < DBL_MAX)
			x2 *= x_max;
	}
	return result;
}


template<int N>
bool Polynomial<N>::operator == (const Polynomial& p) const {
	for (int i = 0; i <= N; i++) { 
		if (coeff_[i] != p.coefficient(i))
			return false; 
	}
	return true;
}


template<int N>
bool Polynomial<N>::operator != (const Polynomial& p) const {
	for (int i = 0; i <= N; i++) { 
		if (coeff_[i] != p.coefficient(i))
			return true; 
	}
	return false;
}


template<int N>
bool Polynomial<N>::is_zero(void) const {
	for (int i = 0; i <= N; i++) { 
		if (coeff_[i] != 0)
			return false; 
	}
	return true;
}


template<int N>
void Polynomial<N>::set_zero(void) {
	memset(coeff_, 0, sizeof(double) * nb_coeffs);
}

template<int N>
Polynomial<N>& Polynomial<N>::add_scaled(const Polynomial& p, double s){
	for (int i = 0; i <= N; i++) {
		coeff_[i] += p.coefficient(i) * s;
	}
	return *this;
}


template<int N>
Polynomial<N>& Polynomial<N>::operator += (const Polynomial<N>& p){
	for (int i = 0; i <= N; i++) {
		coeff_[i] += p.coefficient(i);
	}
	return *this;
}


template<int N>
Polynomial<N>& Polynomial<N>::operator -= (const Polynomial<N>& p) {
	for (int i = 0; i <= N; i++) {
		coeff_[i] -= p.coefficient(i);
	}
	return *this;
}


template<int N>
Polynomial<N> Polynomial<N>::operator + (const Polynomial<N>& p) const{
	Polynomial q;
	for (int i = 0; i <= N; i++) {
		q.coefficient(i) = (coeff_[i] + p.coefficient(i));
	}
	return q;
}


template<int N>
Polynomial<N> Polynomial<N>::operator - (const Polynomial<N>& p) const{
	Polynomial q;
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = (coeff_[i] - p.coefficient(i)); 
	}
	return q;
}


template<int N>
void Polynomial<N>::scale(const Polynomial& p, double w, Polynomial& q){
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = p.coefficient(i) * w; 
	}
}


template<int N>
void Polynomial<N>::add_scaled(const Polynomial& p1, double w1, const Polynomial& p2, double w2, Polynomial& q){
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = p1.coefficient(i) * w1 + p2.coefficient(i) * w2; 
	}
}


template<int N>
void Polynomial<N>::add_scaled(const Polynomial& p1, double w1, const Polynomial& p2, Polynomial& q){
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = p1.coefficient(i) * w1 + p2.coefficient(i); 
	}
}


template<int N>
void Polynomial<N>::add_scaled(const Polynomial& p1, const Polynomial& p2, double w2, Polynomial& q){
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = p1.coefficient(i) + p2.coefficient(i) * w2; 
	}
}


template<int N>
void Polynomial<N>::subtract(const Polynomial &p1, const Polynomial& p2, Polynomial& q){
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = p1.coefficient(i) - p2.coefficient(i); 
	}
}


template<int N>
void Polynomial<N>::negate(const Polynomial& in, Polynomial& out){
	out = in;
	for (int i = 0; i <= N; i++) { 
		out.coefficient(i) = -out.coefficient(i); 
	}
}


template<int N>
Polynomial<N> Polynomial<N>::operator - (void) const{
	Polynomial q = *this;
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = -q.coefficient(i); 
	}
	return q;
}


template<int N>
template<int N2>
Polynomial<N + N2> Polynomial<N>::operator * (const Polynomial<N2>& p) const{
	Polynomial<N + N2> q;
	for (int i = 0; i <= N; i++) { 
		for (int j = 0; j <= N2; j++) { 
			q.coefficient(i + j) += coeff_[i] * p.coefficient(j); 
		} 
	}
	return q;
}

template<int N>
Polynomial<N>& Polynomial<N>::operator += (double s)
{
	coeff_[0] += s;
	return *this;
}


template<int N>
Polynomial<N>& Polynomial<N>::operator -= (double s)
{
	coeff_[0] -= s;
	return *this;
}


template<int N>
Polynomial<N>& Polynomial<N>::operator *= (double s)
{
	for (int i = 0; i <= N; i++) { 
		coeff_[i] *= s; 
	}
	return *this;
}


template<int N>
Polynomial<N>& Polynomial<N>::operator /= (double s)
{
	for (int i = 0; i <= N; i++) { 
		coeff_[i] /= s;
	}
	return *this;
}


template<int N>
Polynomial<N> Polynomial<N>::operator + (double s) const
{
	Polynomial<N> q = *this;
	q.coefficient(0) += s;
	return q;
}


template<int N>
Polynomial<N> Polynomial<N>::operator - (double s) const
{
	Polynomial q = *this;
	q.coefficient(0) -= s;
	return q;
}


template<int N>
Polynomial<N> Polynomial<N>::operator * (double s) const
{
	Polynomial q;
	for (int i = 0; i <= N; i++) { 
		q.coefficient(i) = coeff_[i] * s; 
	}
	return q;
}


template<int N>
Polynomial<N> Polynomial<N>::operator / (double s) const
{
	Polynomial q;
	for (int i = 0; i <= N; i++) {
		q.coefficient(i) = coeff_[i] / s;
	}
	return q;
}


template<int N>
Polynomial<N> Polynomial<N>::scale(double s) const
{
	Polynomial q = *this;
	double s2 = 1.0;
	for (int i = 0; i <= N; i++){
		q.coefficient(i) *= s2;
		s2 /= s;
	}
	return q;
}


template<int N>
Polynomial<N> Polynomial<N>::shift(double t) const
{
	Polynomial<N> q;
	for (int i = 0; i <= N; i++){
		double temp = 1;
		for (int j = i; j >= 0; j--){
			q.coefficient(j) += coeff_[i] * temp;
			temp *= -t*j;
			temp /= (i - j + 1);
		}
	}
	return q;
}


template<int N>
void Polynomial<N>::print(void) const{
	for (int j = 0; j <= N; j++) {
		printf("%6.4f x^%d", coefficients()[j], j);
		if (j < N && coefficients()[j + 1] >= 0) {
			printf(" + ");
		}
	}
	printf("\n");
}


template<int N>
void Polynomial<N>::get_solutions(double value, std::vector<double>& roots, double eps) const
{
	double r[4][2];
	int rCount = 0;
	roots.clear();
	switch (N){
	case 1:
		rCount = Factor(coeff_[1], coeff_[0] - value, r, eps);
		break;
	case 2:
		rCount = Factor(coeff_[2], coeff_[1], coeff_[0] - value, r, eps);
		break;
	case 3:
		rCount = Factor(coeff_[3], coeff_[2], coeff_[1], coeff_[0] - value, r, eps);
		break;
	case 4:
		rCount = Factor(coeff_[4], coeff_[3], coeff_[2], coeff_[1], coeff_[0] - value, r, eps);
		break;
	default:
		printf("Can't solve polynomial of degree: %d\n", N);
	}
	for (int i = 0; i < rCount; i++) {
		if (fabs(r[i][1]) <= eps) {
			roots.push_back(r[i][0]);
		}
	}
}


template< int N >
int Polynomial<N>::get_solutions(double value, double* roots, double eps) const
{
	double _roots[4][2];
	int    _rCount = 0;
	switch (N)
	{
	case 1: 
		_rCount = Factor(coeff_[1], coeff_[0] - value, _roots, eps);
		break;
	case 2: 
		_rCount = Factor(coeff_[2], coeff_[1], coeff_[0] - value, _roots, eps);
		break;
	case 3: 
		_rCount = Factor(coeff_[3], coeff_[2], coeff_[1], coeff_[0] - value, _roots, eps);
		break;
	case 4: 
		_rCount = Factor(coeff_[4], coeff_[3], coeff_[2], coeff_[1], coeff_[0] - value, _roots, eps);
		break;
	default: 
		printf("Can't solve polynomial of degree: %d\n", N);
	}
	int rCount = 0;
	for (int i = 0; i < _rCount; i++) {
		if (fabs(_roots[i][1]) <= eps)
			roots[rCount++] = _roots[i][0];
	}
	return rCount;
}


template< >
Polynomial< 0 > Polynomial< 0 >::BSpline_component(int i)
{
	Polynomial p;
	p.coeff_[0] = 1.;
	return p;
}


template< int N >
Polynomial< N > Polynomial< N >::BSpline_component(int i)
{
	Polynomial p;
	if (i > 0)
	{
		Polynomial< N > _p = Polynomial< N - 1 >::BSpline_component(i - 1).integral();
		p -= _p;
		p.coefficient(0) += _p(1);
	}
	if (i < N)
	{
		Polynomial< N > _p = Polynomial< N - 1 >::BSpline_component(i).integral();
		p += _p;
	}
	return p;
}



#endif


