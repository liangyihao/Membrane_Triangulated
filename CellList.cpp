#include "CellList.hpp"
CellList::CellList(double L,double radius) {
	this->L=L;
	this->radius=radius;
	dL=2*radius;
	NC=L/dL+2;
	X.resize(0);
	InWhichCell.resize(0);
	Cells.resize(0);

	//now init Cells
	Cell emptyCell;
	for(int k=0;k<NC*NC*NC;k++)Cells.push_back(emptyCell);

	//now load coordinates
	double3 Xtemp;
	ifstream in("Pos.txt");
	int3 temp;
	temp.y=-1;
	while(!in.eof()) {
		in>>temp.x>>Xtemp.x>>Xtemp.y>>Xtemp.z;
		if(temp.y==temp.x)break;
		temp.y=temp.x;
		X.push_back(Xtemp);
	}
	in.close();

	//Renormalize length scale
	double mindist=L,maxdist=0;
	double3 dX;
	double dr;
	ifstream In("Tri.txt");
	int f,v1,v2,v3;
	int oldF=-1;
	while(!In.eof()) {
		In>>f>>v1>>v2>>v3;
		if(oldF==f)break;//To prevent the re-reading at the end of file
		dX=X[v1]-X[v2];
		dr=sqrt(DotProd(dX,dX));
		if(dr<mindist)mindist=dr;if(dr>maxdist)maxdist=dr;

		dX=X[v2]-X[v3];
		dr=sqrt(DotProd(dX,dX));
		if(dr<mindist)mindist=dr;if(dr>maxdist)maxdist=dr;

		dX=X[v3]-X[v1];
		dr=sqrt(DotProd(dX,dX));
		if(dr<mindist)mindist=dr;if(dr>maxdist)maxdist=dr;
	}
	In.close();
	if(radius/mindist>sqrt(3)*radius/maxdist){cout<<"Wrong constructing"<<endl;exit(0);}
	double A;
	A=radius/mindist + sqrt(3)*radius/maxdist;
	for(int I=0;I<X.size();I++) X[I]=X[I]/(1.0/A);

	///////////////////////////////////Make Cell List///////////////////
	for(int I=0;I<X.size();I++) {
		Xtemp=X[I];
		temp.x=(Xtemp.x+L/2)/dL;
		temp.y=(Xtemp.y+L/2)/dL;
		temp.z=(Xtemp.z+L/2)/dL;
		Cells[temp.z*NC*NC+temp.y*NC+temp.x].Insert(I);
		InWhichCell.push_back(temp);
	}
}

bool CellList::Overlap(int I,double3 NX) {
	if(NX.x<-L/2)return true;
	if(NX.x>+L/2)return true;

	if(NX.y<-L/2)return true;
	if(NX.y>+L/2)return true;

	if(NX.z<-L/2)return true;
	if(NX.z>+L/2)return true;

	int i0,j0,k0;
	i0=(NX.x+L/2)/dL;
	j0=(NX.y+L/2)/dL;
	k0=(NX.z+L/2)/dL;
	//double dx,dy,dz,dr2;
	double3 dX;
	int J;
	for(int i=max(i0-1,0);i<=min(i0+1,NC-1);i++)
		for(int j=max(j0-1,0);j<=min(j0+1,NC-1);j++)
			for(int k=max(k0-1,0);k<=min(k0+1,NC-1);k++)
				for(int l=0;l<Cells[k*NC*NC+j*NC+i].num_particle;l++) {
					J=Cells[k*NC*NC+j*NC+i].particle[l];
					if(I==J)continue;
					dX=X[J]-NX;
					if(DotProd(dX,dX)<4*radius*radius)return true;
				}
	return false;
}

void CellList::Update(int I,double3 NX) {
	if((NX.x<-L/2)||(NX.x>+L/2)){cout<<"Error in updating, new position out of boundary"<<endl;return;}
	if((NX.y<-L/2)||(NX.y>+L/2)){cout<<"Error in updating, new position out of boundary"<<endl;return;}
	if((NX.z<-L/2)||(NX.z>+L/2)){cout<<"Error in updating, new position out of boundary"<<endl;return;}
	int3 tempInt3;
	tempInt3=InWhichCell[I];
	Cells[tempInt3.z*NC*NC+tempInt3.y*NC+tempInt3.x].Delete(I);

	X[I]=NX;
	tempInt3.x=(NX.x+L/2)/dL;
	tempInt3.y=(NX.y+L/2)/dL;
	tempInt3.z=(NX.z+L/2)/dL;
	InWhichCell[I]=tempInt3;
	Cells[tempInt3.z*NC*NC+tempInt3.y*NC+tempInt3.x].Insert(I);
}
