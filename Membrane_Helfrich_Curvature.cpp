#include "Membrane_Helfrich.hpp"
double Membrane_Helfrich::Mean_Curvature_Itzykson(int v){//Should be checked
	/*Based on the paper:
	ItzyksonC1986Proc.oftheGIFTSeminar,Jaca85edJAbad et al pp 130–88
	
	Gaussian Curvature:The Gauss-Bonnet Scheme
 	*/
 	double3 H;
 	H.x=0;H.y=0;H.z=0;
 	double lambda=0;
 	double temp;

	double cost1,cost2,ctgt1,ctgt2;
 	double3 dx1,dx2,dx3,dr1,dr2;
 	double dx1l,dx2l,dx3l,dr1l,dr2l;
 	int num_nei;
 	num_nei=Connection[v].size();
 	dx1=X[Connection[v][0].x]-X[v];
 	dx2=X[Connection[v][num_nei-1].x]-X[v];
 	dx3=X[Connection[v][num_nei-2].x]-X[v];
 	dr1=dx1-dx2;
 	dr2=dx2-dx3;
 	dx1l=sqrt(DotProd(dx1,dx1));
 	dx2l=sqrt(DotProd(dx2,dx2));
 	dx3l=sqrt(DotProd(dx3,dx3));
 	dr1l=sqrt(DotProd(dr1,dr1));
 	dr2l=sqrt(DotProd(dr2,dr2));

 	//cout<<"Computing Mean Curvature on the node "<<v<<":"<<endl<<"Number of neighbors: "<<num_nei<<endl;
 	for(int k=0;k<Connection[v].size();k++) {
 		//cout<<"neighbor "<<k<<endl;
 		//cout<<"dx1: "<<dx1.x<<' '<<dx1.y<<' '<<dx1.z<<' '<<dx1l<<endl;
 		//cout<<"dx2: "<<dx2.x<<' '<<dx2.y<<' '<<dx2.z<<' '<<dx2l<<endl;
 		//cout<<"dx3: "<<dx3.x<<' '<<dx3.y<<' '<<dx3.z<<' '<<dx3l<<endl;
 		//cout<<"dr1: "<<dr1.x<<' '<<dr1.y<<' '<<dr1.z<<' '<<dr1l<<endl;
 		//cout<<"dr2: "<<dr2.x<<' '<<dr2.y<<' '<<dr2.z<<' '<<dr2l<<endl;

 		cost1=DotProd(dr1,dx1)/(dr1l*dx1l);
 		ctgt1=cost1/sqrt(1-cost1*cost1);
 		cost2=-DotProd(dr2,dx3)/(dr2l*dx3l);
 		ctgt2=cost2/sqrt(1-cost2*cost2);

 		temp=(ctgt1+ctgt2)/2.0;
 		H = H + temp*dx2;
 		lambda +=0.25*temp*dx2l*dx2l;

 		dx3=dx2;
 		dx2=dx1;
 		dx1=X[Connection[v][(k+1)%num_nei].x]-X[v];
 		dr2=dr1;
 		dr1=dx1-dx2;
 		dx3l=dx2l;
 		dx2l=dx1l;
	 	dx1l=sqrt(DotProd(dx1,dx1));
	 	dr2l=dr1l;
 		dr1l=sqrt(DotProd(dr1,dr1));

 	}

 	H = H/(2*lambda);
 	int sign;
 	if(DotProd(H,Unit_Normal_of_Triangle[Faces_of_Vertex[v][0]])>0)sign=-1;else sign=1;
 	return sign*sqrt(DotProd(H,H));
}


double Membrane_Helfrich::Gauss_Curvature_STD(int v){//Should be checked
 	//Step4: compute Gaussian Curvature on v
	double sum_alpha=0,sum_A=0;
	double3 dx,dx_old;
	dx=X[v]-X[Connection[v].back().x];
	for(int k=0;k<Connection[v].size();k++) {
		dx_old=dx;
		dx=X[v]-X[Connection[v][k].x];
		sum_alpha+=acos(DotProd(dx,dx_old)/sqrt(DotProd(dx,dx)*DotProd(dx_old,dx_old)));
	}

	for(int k=0;k<Faces_of_Vertex[v].size();k++) {
		int f;
		f=Faces_of_Vertex[v][k];
		sum_A+=Area_of_Triangle[f];
	}
	return 3*(2*M_PI-sum_alpha)/sum_A;	
}

