
#ifndef _MATH_FACTOR_H_
#define _MATH_FACTOR_H_

#include "math_common.h"


#define PI		3.1415926535897932384
#define SQRT_3	1.7320508075688772935

double MATH_API ArcTan2(double y,double x);
double MATH_API Angle(const double in[2]);
void MATH_API Sqrt(const double in[2], double out[2]);
void MATH_API Add(const double in1[2], const double in2[2], double out[2]);
void MATH_API Subtract(const double in1[2], const double in2[2], double out[2]);
void MATH_API Multiply(const double in1[2], const double in2[2], double out[2]);
void MATH_API Divide(const double in1[2], const double in2[2], double out[2]);

int MATH_API Factor(double a1, double a0, double roots[1][2], double EPS);
int MATH_API Factor(double a2, double a1, double a0, double roots[2][2], double EPS);
int MATH_API Factor(double a3, double a2, double a1, double a0, double roots[3][2], double EPS);
int MATH_API Factor(double a4, double a3, double a2, double a1, double a0, double roots[4][2], double EPS);

int MATH_API Solve(const double* eqns, const double* values, double* solutions, int dim);

#endif // _MATH_FACTOR_H_