#ifndef PUBLIC_H
#define PUBLIC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;
typedef struct {
	int x,y;
}int2;
typedef struct {
	int x,y,z;
}int3;

typedef struct {
	double x,y,z;
}double3;

double3 operator+(double3 const&A1,double3 const&A2);
double3 operator-(double3 const&A1,double3 const&A2);
double3 operator/(double3 const&A,double const&lambda);
double3 operator*(double const&lambda,double3 const&A);
double3 CrossProd(double3 const&A1,double3 const&A2);//cross product
double DotProd(double3 const&A1,double3 const&A2);//dot product

#define MINDEGREE 5
#define MAXDEGREE 7
#endif
