#include "TriMesh.hpp"
#include "CellList.hpp"
#include "Rand4.hpp"
class Membrane_Helfrich: public TriMesh, public CellList {
 private:
	double K,H_0,K_G;
	//double Ka,A0;
	//double dP;//Osmotic Pressure
 	double Ka,Surface_A0,Surface_A;//These parameters are used to control the total area with Energy 0.5*Ka*(1-A/A0)^2
 	double Kv,Cell_V0,Cell_V;//These parameters are for controling the total volumn with Energy 0.5*Kv*(1-V/V0)^2

 	double MOVER;
 	vector<double3>Center_of_Triangle;
 	vector<double>Area_of_Triangle;
 	vector<double3>Unit_Normal_of_Triangle;
 	vector<double>Determinant_of_Triangle_over6;//For the computation of Volumn

 	vector<double>Mean_Curvature;
 	vector<double>Gauss_Curvature;
 	void Compute_Center_of_Triangle(int f);
 	void Compute_Area_Unit_Normal_of_Triangle(int f);
 	void Compute_Det_of_Triangle_over6(int f);//For the computation of Volumn
	
	double Mean_Curvature_Itzykson(int v);//Should be checked
	double Gauss_Curvature_STD(int v);//Should be checked
	void UpdateGeoInfo(int v) {
		 	Mean_Curvature[v]=Mean_Curvature_Itzykson(v);
 			Gauss_Curvature[v]=Gauss_Curvature_STD(v);
	}

	double Helfrich_Energy(int v);

 public:
	Membrane_Helfrich(double L,double radius,double MOVER,double K,double H_0,double K_G,double Ka,double Kv,double Surface_A0, double eta):CellList(L,radius),MOVER(MOVER),K(K),H_0(H_0),K_G(K_G),Ka(Ka),Kv(Kv),Surface_A0(Surface_A0) {
 		Center_of_Triangle.resize(F);
 		Area_of_Triangle.resize(F);
 		Unit_Normal_of_Triangle.resize(F);
 		Determinant_of_Triangle_over6.resize(F);
 		Mean_Curvature.resize(V);
 		Gauss_Curvature.resize(V);

 		for(int f=0;f<F;f++) {
 			Compute_Center_of_Triangle(f);
 			Compute_Area_Unit_Normal_of_Triangle(f);
 			Compute_Det_of_Triangle_over6(f);
 		}

 		for(int v=0;v<V;v++)UpdateGeoInfo(v);
		srand4(0);

		
		//Check dependence.
		if(eta>1){cout<<"Error eta should not be greater than 1"<<endl;exit(0);}
		if(1.0/abs(H_0)<10*radius){cout<<"Triangle too large, try to increase number of triangles"<<endl;exit(0);}
		
		if(Surface_A0<F*radius*radius*sqrt(3.0)){cout<<"Total area too small, the current settings cannot approach. Total are should between "<< F*radius*radius*sqrt(3.0)<<" and "<< 3*F*radius*radius*sqrt(3.0)<<endl;exit(0);}
		if(Surface_A0>3*F*radius*radius*sqrt(3.0)){cout<<"Total area too large, the current settings cannot approach. Total are should between "<< F*radius*radius*sqrt(3.0)<<" and "<< 3*F*radius*radius*sqrt(3.0)<<endl;exit(0);}
				

		//Check environment
		cout<<"Variance of Area is about "<<Surface_A0/sqrt(Ka)<<endl;
		cout<<"Variance of Volume is about "<<Cell_V0*sqrt(2/Kv)<<endl;
		cout<<"Variance of mean curve is about "<<1.0/sqrt(K)<<endl;
		//if(Ka<K*F*F)cout<<"Warning Ka too small, constraint of surface area may fail"<<endl;
		//if(Ka<K_G*F*F)cout<<"Warning Ka too small, constraint of surface area may fail"<<endl;
		//if(Kv<K*F*F)cout<<"Warning K_G too small, constraint of volume may fail"<<endl;
		//if(Kv<K_G*F*F)cout<<"Warning K_G too small, constraint of volume may fail"<<endl;

		//Set size para
		Cell_V0=eta*4*M_PI/3*pow(Surface_A0/(4*M_PI),1.5);
		Cell_V=0;
		Surface_A=0;
		for(int f=0;f<F;f++) {
			Surface_A+=Area_of_Triangle[f];
			Cell_V+=Determinant_of_Triangle_over6[f];
		}
 	}

 	bool move_node(int I);//true:accepted, false:rejected
 	bool flip_edge(int e);
	

	void Print(int filenum);
	void Check();
	
 	void SetTrialPara(double MOVER){
 		cout<<"Set MOVER:"<<(this->MOVER)<<"->"<<MOVER<<endl;
 		this->MOVER=MOVER;
 	}

	void SetMembranePara(double K,double K_G,double Ka,double Kv) {
 		this->K=K;
 		this->K_G=K_G;
 		this->Ka=Ka;
 		this->Kv=Kv;
		//Check environment
		cout<<"Variance of Area is about "<<Surface_A0/sqrt(Ka)<<endl;
		cout<<"Variance of Volume is about "<<Cell_V0*sqrt(2/Kv)<<endl;
		cout<<"Variance of mean curve is about "<<1.0/sqrt(K)<<endl;

		//if(Ka/K<F*F)cout<<"Warning Ka too small, constraint of surface area may fail"<<endl;
		//if(Ka/K_G<F*F)cout<<"Warning Ka too small, constraint of surface area may fail"<<endl;
		//if(Kv/K<F*F)cout<<"Warning K_G too small, constraint of volume may fail"<<endl;
		//if(Kv/K_G<F*F)cout<<"Warning K_G too small, constraint of volume may fail"<<endl;

 	}

 	void Refresh_Sys_Var() {
		Cell_V=0;
		Surface_A=0;
		for(int f=0;f<F;f++) {
			Surface_A+=Area_of_Triangle[f];
			Cell_V+=Determinant_of_Triangle_over6[f];
		} 		
 	}
};
