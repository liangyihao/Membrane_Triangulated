#include "Membrane_Helfrich.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>

void Membrane_Helfrich::Compute_Center_of_Triangle(int f) {
	double3 X1,X2,X3;
	int v1,v2,v3;
	v1=Vertices_of_Face[f][0];
	v2=Vertices_of_Face[f][1];
	v3=Vertices_of_Face[f][2];
	X1=X[v1];X2=X[v2];X3=X[v3];
	Center_of_Triangle[f]=(X1+X2+X3)/3;
}

void Membrane_Helfrich::Compute_Area_Unit_Normal_of_Triangle(int f) {
	double3 X1,X2,X3,Xc;
	int v1,v2,v3;
	double3 dr1,dr2;
	v1=Vertices_of_Face[f][0];
	v2=Vertices_of_Face[f][1];
	v3=Vertices_of_Face[f][2];
	X1=X[v1];X2=X[v2];X3=X[v3];
	dr1=X2-X1;
	dr2=X3-X1;
	
	double3 N;
	double A;
	N=CrossProd(dr1,dr2);
	A=sqrt(DotProd(N,N));
	Area_of_Triangle[f]=0.5*A;
	Unit_Normal_of_Triangle[f]=N/A;
}

void Membrane_Helfrich::Compute_Det_of_Triangle_over6(int f) {//For the computation of Volumn
	double3 r1,r2,r3;
	r1=X[Vertices_of_Face[f][0]];
	r2=X[Vertices_of_Face[f][1]];
	r3=X[Vertices_of_Face[f][2]];
	Determinant_of_Triangle_over6[f]=(r1.x*(r2.y*r3.z-r2.z*r3.y)   - r1.y*(r2.x*r3.z-r2.z*r3.x)    + r1.z*(r2.x*r3.y-r2.y*r3.x))/6.0;
}

double Membrane_Helfrich::Helfrich_Energy(int v){
	double A=0;
	int f;
	for(int k=0;k<Faces_of_Vertex[v].size();k++){
		f=Faces_of_Vertex[v][k];
		A+=Area_of_Triangle[f];
	}
	A/=3.0;
	double ans=0,temp;
	temp=Mean_Curvature[v]-H_0;
	ans=0.5*K*temp*temp+K_G*Gauss_Curvature[v];
	return ans*A;
}

void Membrane_Helfrich::Check() {//Rewrite
	cout<<"Check"<<endl;
	cout<<"Area="<<Surface_A<<" Std Area="<<Surface_A0<<endl;
	cout<<"Volume="<<Cell_V<<" Std Volume="<<Cell_V0<<endl;
	double3 T1,T2;
	double temp1,temp2;
	for(int f=0;f<F;f++) {
		T1=Center_of_Triangle[f];
		T2=Unit_Normal_of_Triangle[f];
		temp1=Area_of_Triangle[f];
		temp2=Determinant_of_Triangle_over6[f];
		Compute_Center_of_Triangle(f);
		Compute_Area_Unit_Normal_of_Triangle(f);
		Compute_Det_of_Triangle_over6(f);

		T1 = T1 - Center_of_Triangle[f];
		T2 = T2 -Unit_Normal_of_Triangle[f];
		temp1 = temp1- Area_of_Triangle[f];
		temp2 = temp2 -Determinant_of_Triangle_over6[f];
		if(DotProd(T1,T1)>(1E-10))cout<<"Warning Result does not match"<<endl;
		if(DotProd(T2,T2)>(1E-10))cout<<"Warning Result does not match"<<endl;
		if(temp1*temp1>(1E-10))cout<<"Warning Result does not match"<<endl;
		if(temp2*temp2>(1E-10))cout<<"Warning Result does not match"<<endl;
 	}
	for(int v=0;v<V;v++) {
		temp1=Mean_Curvature[v]-Mean_Curvature_Itzykson(v);
		temp2=Gauss_Curvature[v]-Gauss_Curvature_STD(v);

		if(temp1*temp1>(1E-10))cout<<"Warning Result does not match"<<endl;
		if(temp2*temp2>(1E-10))cout<<"Warning Result does not match"<<endl;
	}

	double tempA=0,tempV=0;
	for(int f=0;f<F;f++) {
		tempA+=Area_of_Triangle[f];
		tempV+=Determinant_of_Triangle_over6[f];
	}
	if(abs(tempA-Surface_A)/F>1E-10)cout<<"Warning Result does not match in area "<<tempA<<" "<<Surface_A<<endl;
	if(abs(tempV-Cell_V)/F>1E-10)cout<<"Warning Result does not match in Volume "<<tempV<<" "<<Cell_V<<endl;
}


void Membrane_Helfrich::Print(int filenum) {
	string TXT(".txt");
	string FNUM;
	FNUM=to_string(filenum);

	string FF("Tri");
	FF+=FNUM;
	FF+=TXT;
	ofstream Fout(FF);
	for(int k=0;k<F;k++) {
		Fout<<k;
		for(int l=0;l<Vertices_of_Face[k].size();l++)Fout<<' '<<Vertices_of_Face[k][l];
		Fout<<endl;
	}
	Fout.close();

	string VF("Pos");
	VF+=FNUM;
	VF+=TXT;
	ofstream Vout(VF);
	for(int k=0;k<V;k++) Vout<<k<<' '<<X[k].x<<' '<<X[k].y<<' '<<X[k].z<<endl;
	Vout.close();
	
	double MeanC,GaussC,C1,C2;
	double a,b,c,SqrtDelta;
	string CV("Curvature");
	CV+=FNUM;
	CV+=TXT;
	ofstream CVout(CV);
	for(int v=0;v<V;v++) {
		MeanC=Mean_Curvature[v];
		GaussC=Gauss_Curvature[v];
		a=1;
		b=2*MeanC;
		c=GaussC;
		SqrtDelta=sqrt(b*b-4*a*c);
		C1=(-b+SqrtDelta)/2;
		C2=(-b-SqrtDelta)/2;
		//CVout<<v<<' '<<MeanC<<' '<<GaussC<<' '<<1.0/C1<<' '<<1.0/C2<<endl;
		CVout<<v<<' '<<MeanC<<' '<<GaussC<<endl;
	}
	CVout.close();
	//cout<<""

	double Etot=0,temp;
	for(int k=0;k<V;k++)Etot+=Helfrich_Energy(k);
	cout<<"Helfrich Energy: "<<Etot<<endl;

	Etot+=0.5*Ka*pow(1.0-Surface_A/Surface_A0,2.0);
	Etot+=0.5*Kv*pow(1.0-Cell_V/Cell_V0,2.0);
	cout<<"Surface tension energy "<<0.5*Ka*pow(1.0-Surface_A/Surface_A0,2.0)<<endl;
	cout<<"Volume energy "<<0.5*Kv*pow(1.0-Cell_V/Cell_V0,2.0)<<endl;
	cout<<"Overall Energy: "<<Etot<<endl;fflush(stdout);

}
