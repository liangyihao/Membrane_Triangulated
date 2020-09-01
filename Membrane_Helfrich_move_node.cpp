#include "Membrane_Helfrich.hpp"
bool Membrane_Helfrich::move_node(int I) {
	double3 dr;
	int num_nei=Faces_of_Vertex[I].size();
	//////////////////Propose a new position and check whether it overlaps with others or longger than the tether length
	double3 dX,NewX,OldX;
	dX.x=(rand4()-0.5)*MOVER;dX.y=(rand4()-0.5)*MOVER;dX.z=(rand4()-0.5)*MOVER;
	OldX=X[I];
	NewX=X[I]+dX;
	if(Overlap(I,NewX))return false;
	for(int k=0;k<num_nei;k++) {
		dr=NewX-X[Connection[I][k].x];
		if((DotProd(dr,dr)>12*radius*radius))return false;
	}

	//////////////////Prepare the energy computation//////////////////////////////////////////////////////////////////////

	//////////////////Compute Old Energy//////////////////////////////////////////////////////////////////////////////////
	double E_old=0;
	E_old+=Helfrich_Energy(I);
	for(int k=0;k<num_nei;k++) {
		E_old+=Helfrich_Energy(Connection[I][k].x);
	}

	E_old+=0.5*Ka*pow(1.0-Surface_A/Surface_A0,2.0);
	E_old+=0.5*Kv*pow(1.0-Cell_V/Cell_V0,2.0);

	///////////////////!!!!!!!!!!!!!!!!!!///////////////////////////////////////////////////////////////////////////////

	//////////////////Now propose a new position and update new geometrical values///////////////////////////////////
	double TempArea1=0,TempArea2=0,TempV1=0,TempV2=0;

	Update(I,NewX);
	for(int k=0;k<num_nei;k++) {
		
		TempArea1+=Area_of_Triangle[Faces_of_Vertex[I][k]];
		TempV1+=Determinant_of_Triangle_over6[Faces_of_Vertex[I][k]];

		Compute_Center_of_Triangle(Faces_of_Vertex[I][k]);
		Compute_Area_Unit_Normal_of_Triangle(Faces_of_Vertex[I][k]);
		Compute_Det_of_Triangle_over6(Faces_of_Vertex[I][k]);
		
		TempArea2+=Area_of_Triangle[Faces_of_Vertex[I][k]];
		TempV2+=Determinant_of_Triangle_over6[Faces_of_Vertex[I][k]];
	}
	for(int k=0;k<num_nei;k++) {
		UpdateGeoInfo(I);
		UpdateGeoInfo(Connection[I][k].x);
	}

	//////////////////Compute New Energy/////////////////////////////////////////////////////////////////////////////
	double E_new=0;
	E_new+=Helfrich_Energy(I);
	for(int k=0;k<num_nei;k++) {
		E_new+=Helfrich_Energy(Connection[I][k].x);
	}

	E_new+=0.5*Ka*pow(1.0-(Surface_A + TempArea2 - TempArea1)/Surface_A0,2.0);
	E_new+=0.5*Kv*pow(1.0-(Cell_V + TempV2 - TempV1)/Cell_V0,2.0);


	//////////////////Make decision//////////////////////////////////////////////////////////////////////////////////
	if(rand4()>exp(E_old-E_new)){
		Update(I,OldX);
		for(int k=0;k<num_nei;k++) {
			Compute_Center_of_Triangle(Faces_of_Vertex[I][k]);
			Compute_Area_Unit_Normal_of_Triangle(Faces_of_Vertex[I][k]);
			Compute_Det_of_Triangle_over6(Faces_of_Vertex[I][k]);
		}
		for(int k=0;k<num_nei;k++) {
			UpdateGeoInfo(I);
			UpdateGeoInfo(Connection[I][k].x);
		}
		return false;
	}

	Surface_A+=(TempArea2-TempArea1);
	Cell_V+=(TempV2-TempV1);
	return true;
}
