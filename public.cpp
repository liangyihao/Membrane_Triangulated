#include "public.hpp"
double3 operator+(double3 const&A1,double3 const&A2) {
	double3 A;
	A.x=A1.x+A2.x;
	A.y=A1.y+A2.y;
	A.z=A1.z+A2.z;
	return A;
}

double3 operator-(double3 const&A1,double3 const&A2) {
	double3 A;
	A.x=A1.x-A2.x;
	A.y=A1.y-A2.y;
	A.z=A1.z-A2.z;
	return A;
}

double3 operator/(double3 const&A,double const&lambda) {
	double3 As;
	As.x=A.x/lambda;
	As.y=A.y/lambda;
	As.z=A.z/lambda;
	return As;
}

double3 operator*(double const&lambda,double3 const&A) {
	double3 As;
	As.x=A.x*lambda;
	As.y=A.y*lambda;
	As.z=A.z*lambda;
	return As;
}

double3 CrossProd(double3 const&A1,double3 const&A2) {//cross product
	double3 A;
	A.x=A1.y*A2.z-A1.z*A2.y;
	A.y=A1.z*A2.x-A1.x*A2.z;
	A.z=A1.x*A2.y-A1.y*A2.x;
	return A;
}

double DotProd(double3 const&A1,double3 const&A2){//dot product
	return A1.x*A2.x + A1.y*A2.y + A1.z*A2.z;
}
