/*此段需要严格进行测试，注意拓扑关系*/
#include "Membrane_Helfrich.hpp"
bool Membrane_Helfrich::flip_edge(int e) {
	////pre-reject
	int v1,v2,u1,u2;
	double3 dr;
	v1=Vertices_of_Edge[e][0];
	v2=Vertices_of_Edge[e][1];
	u1=ExVertices_of_AdjFaces_of_Edge[e][0];
	u2=ExVertices_of_AdjFaces_of_Edge[e][1];
	if(Connection[v1].size()==MINDEGREE)return false;
	if(Connection[v2].size()==MINDEGREE)return false;
	if(Connection[u1].size()==MAXDEGREE)return false;
	if(Connection[u2].size()==MAXDEGREE)return false;

	dr=X[u1]-X[u2];
	if(DotProd(dr,dr)>12*radius*radius)return false;


	int f1,f2;//注意!未接受的话需要再进行一次flip,其结果是f1,f2编号顺序会颠倒
	f1=AdjFaces_of_Edge[e][0];
	f2=AdjFaces_of_Edge[e][1];
	/////////////////////////////////Compute Energy//////////////////////


	double E_old=0;
	double temp;
	E_old+=Helfrich_Energy(v1);
	E_old+=Helfrich_Energy(v2);
	E_old+=Helfrich_Energy(u1);
	E_old+=Helfrich_Energy(u2);


	E_old+=0.5*Ka*pow(1.0-Surface_A/Surface_A0,2.0);
	E_old+=0.5*Kv*pow(1.0-Cell_V/Cell_V0,2.0);

	///////////////////////////////Flip///////////////////////////////////////应当检查n, N 耦合项不变，或者检查新的cos(theta)
	double TempArea1=0,TempArea2=0,TempV1=0,TempV2=0;
	TempArea1+=Area_of_Triangle[f1];
	TempV1+=Determinant_of_Triangle_over6[f1];
	TempArea1+=Area_of_Triangle[f2];
	TempV1+=Determinant_of_Triangle_over6[f2];

	Flip(e);
	Compute_Center_of_Triangle(f1);
	Compute_Center_of_Triangle(f2);
	Compute_Area_Unit_Normal_of_Triangle(f1);
	Compute_Area_Unit_Normal_of_Triangle(f2);
	Compute_Det_of_Triangle_over6(f1);
	Compute_Det_of_Triangle_over6(f2);
	UpdateGeoInfo(v1);
	UpdateGeoInfo(v2);
	UpdateGeoInfo(u1);
	UpdateGeoInfo(u2);
	
	TempArea2+=Area_of_Triangle[f1];
	TempV2+=Determinant_of_Triangle_over6[f1];
	TempArea2+=Area_of_Triangle[f2];
	TempV2+=Determinant_of_Triangle_over6[f2];


	//////////////////////////Compute New Energy//////////////////////////////
	double E_new=0;	
	E_new+=Helfrich_Energy(v1);
	E_new+=Helfrich_Energy(v2);
	E_new+=Helfrich_Energy(u1);
	E_new+=Helfrich_Energy(u2);


	E_new+=0.5*Ka*pow(1.0-(Surface_A + TempArea2 - TempArea1)/Surface_A0,2.0);
	E_new+=0.5*Kv*pow(1.0-(Cell_V + TempV2 - TempV1)/Cell_V0,2.0);

	////////////////////////////Decision or Recover////////////////////////////
	if(rand4()<exp(E_old-E_new)){
		Surface_A+=(TempArea2-TempArea1);
		Cell_V+=(TempV2-TempV1);
		return true;
	}
	//If rejected, recover
	Flip(e);
	Compute_Center_of_Triangle(f1);
	Compute_Center_of_Triangle(f2);
	Compute_Area_Unit_Normal_of_Triangle(f1);
	Compute_Area_Unit_Normal_of_Triangle(f2);	
	Compute_Det_of_Triangle_over6(f1);
	Compute_Det_of_Triangle_over6(f2);
	UpdateGeoInfo(v1);
	UpdateGeoInfo(v2);
	UpdateGeoInfo(u1);
	UpdateGeoInfo(u2);

	return false;
}
