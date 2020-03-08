
#ifndef _MATH_FACTOR_H_
#define _MATH_FACTOR_H_




#define PI		3.1415926535897932384
#define SQRT_3	1.7320508075688772935

double ArcTan2(double y,double x);
double Angle(const double in[2]);
void Sqrt(const double in[2], double out[2]);
void Add(const double in1[2], const double in2[2], double out[2]);
void Subtract(const double in1[2], const double in2[2], double out[2]);
void Multiply(const double in1[2], const double in2[2], double out[2]);
void Divide(const double in1[2], const double in2[2], double out[2]);

int Factor(double a1, double a0, double roots[1][2], double EPS);
int Factor(double a2, double a1, double a0, double roots[2][2], double EPS);
int Factor(double a3, double a2, double a1, double a0, double roots[3][2], double EPS);
int Factor(double a4, double a3, double a2, double a1, double a0, double roots[4][2], double EPS);

int Solve(const double* eqns, const double* values, double* solutions, int dim);

#endif // _MATH_FACTOR_H_