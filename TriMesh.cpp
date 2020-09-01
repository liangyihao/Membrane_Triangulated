#include "TriMesh.hpp"
vector<int> EmptyIntVec(0);
vector<int2> EmptyInt2Vec(0);
int TriMesh::FindOrCreate_Edge(int i,int j) {//仅可适用于构造函数中
	//If not found, create this edge and update all the relation between edge and vertex
	if(i==j){cout<<"Wrong, encountered loop"<<endl;return -1;}
	int e;
	int l1,l2,m1,m2;
	l1=max(i,j);
	l2=min(i,j);
	for(int k=0;k<Vertices_of_Edge.size();k++) {
		m1=max(Vertices_of_Edge[k][0],Vertices_of_Edge[k][1]);
		m2=min(Vertices_of_Edge[k][0],Vertices_of_Edge[k][1]);
		if((l1==m1)&&(l2==m2))return k;
	}
	
	//Now create new edge (i,j)
	e=Vertices_of_Edge.size();
	Vertices_of_Edge.push_back(EmptyIntVec);
	Vertices_of_Edge[e].push_back(i);
	Vertices_of_Edge[e].push_back(j);

	AdjFaces_of_Edge.push_back(EmptyIntVec);//Just expand the array here, do not record information about faces
	ExVertices_of_AdjFaces_of_Edge.push_back(EmptyIntVec);//Just expand the array here, do not record information about faces

	//Now issue the connection
	while(Connection.size()<=l1) {
		Connection.push_back(EmptyInt2Vec);
		Faces_of_Vertex.push_back(EmptyIntVec);
	}
	int2 temp;
	temp.x=j;temp.y=e;
	Connection[i].push_back(temp);
	temp.y=i;
	Connection[j].push_back(temp);
	return e;
}

TriMesh::TriMesh() {
	//Initialize data structures
	AdjFaces_of_Face.resize(0);
	Edges_of_Face.resize(0);
	Vertices_of_Face.resize(0);

	AdjFaces_of_Edge.resize(0);
	Vertices_of_Edge.resize(0);
	ExVertices_of_AdjFaces_of_Edge.resize(0);

	Connection.resize(0);
	Faces_of_Vertex.resize(0);

	//Read faces
	ifstream In("Tri.txt");
	int f,v1,v2,v3,e;
	int oldF=-1;
	while(!In.eof()) {
		In>>f>>v1>>v2>>v3;//In the input file, the clockwise directions are well defined,i.e. v1->v2->v3->v1...
		if(oldF==f)break;//To prevent the re-reading at the end of file
		if(oldF+1!=f){cout<<"Wrong in construct connection"<<endl;break;}
		oldF=f;
		Vertices_of_Face.push_back(EmptyIntVec);
		Edges_of_Face.push_back(EmptyIntVec);
		Vertices_of_Face[f].push_back(v1);
		Vertices_of_Face[f].push_back(v2);
		Vertices_of_Face[f].push_back(v3);

		e=FindOrCreate_Edge(v1,v2);
		AdjFaces_of_Edge[e].push_back(f);
		ExVertices_of_AdjFaces_of_Edge[e].push_back(v3);
		Edges_of_Face[f].push_back(e);

		e=FindOrCreate_Edge(v2,v3);
		AdjFaces_of_Edge[e].push_back(f);
		ExVertices_of_AdjFaces_of_Edge[e].push_back(v1);
		Edges_of_Face[f].push_back(e);
		
		e=FindOrCreate_Edge(v3,v1);
		AdjFaces_of_Edge[e].push_back(f);
		ExVertices_of_AdjFaces_of_Edge[e].push_back(v2);
		Edges_of_Face[f].push_back(e);

		Faces_of_Vertex[v1].push_back(f);
		Faces_of_Vertex[v2].push_back(f);
		Faces_of_Vertex[v3].push_back(f);//In these three lines, the order of Faces_of_Vertex is not the same as that the system defined

		AdjFaces_of_Face.push_back(EmptyIntVec);
	}
	In.close();
	F=Vertices_of_Face.size();V=Faces_of_Vertex.size();E=AdjFaces_of_Edge.size();
	for(f=0;f<F;f++) {
		e=Edges_of_Face[f][0];
		if(AdjFaces_of_Edge[e][0]!=f)AdjFaces_of_Face[f].push_back(AdjFaces_of_Edge[e][0]);
		if(AdjFaces_of_Edge[e][1]!=f)AdjFaces_of_Face[f].push_back(AdjFaces_of_Edge[e][1]);

		e=Edges_of_Face[f][1];
		if(AdjFaces_of_Edge[e][0]!=f)AdjFaces_of_Face[f].push_back(AdjFaces_of_Edge[e][0]);
		if(AdjFaces_of_Edge[e][1]!=f)AdjFaces_of_Face[f].push_back(AdjFaces_of_Edge[e][1]);

		e=Edges_of_Face[f][2];
		if(AdjFaces_of_Edge[e][0]!=f)AdjFaces_of_Face[f].push_back(AdjFaces_of_Edge[e][0]);
		if(AdjFaces_of_Edge[e][1]!=f)AdjFaces_of_Face[f].push_back(AdjFaces_of_Edge[e][1]);
	}

	//Now make Connection and Faces_of_Vertex in the correct clockwise
	Connection.resize(V);
	int2 tempInt2;
	for(int I=0;I<V;I++) {
		f=Faces_of_Vertex[I][0];
		int l;int num_adjf;
		num_adjf=Faces_of_Vertex[I].size();
		Faces_of_Vertex[I].resize(0);
		for(l=0;l<3;l++)if(Vertices_of_Face[f][l]==I)break;
		v1=Vertices_of_Face[f][(l+1)%3];
		v2=Vertices_of_Face[f][(l+2)%3];
		tempInt2.x=v1;
		tempInt2.y=FindOrCreate_Edge(I,v1);
		Connection[I].resize(0);
		for(int k=0;k<num_adjf;k++) {
			Connection[I].push_back(tempInt2);
			Faces_of_Vertex[I].push_back(f);
			//cout<<"Vertex "<<I<<" "<<tempInt2.x<<' '<<tempInt2.y<<' '<<f<<endl;
			v1=v2;
			e=FindOrCreate_Edge(I,v1);
			if(AdjFaces_of_Edge[e][0]!=f)f=AdjFaces_of_Edge[e][0];else f=AdjFaces_of_Edge[e][1];
			for(l=0;l<3;l++)if(Vertices_of_Face[f][l]==I)break;
			v2=Vertices_of_Face[f][(l+2)%3];
			tempInt2.x=v1;
			tempInt2.y=e;
		}
	}
}

int TriMesh::FindEdge(int i,int j) {
	for(int k=0;k<Connection[i].size();k++)if(Connection[i][k].x==j)return Connection[i][k].y;
	return -1;
}
bool TriMesh::Flip(int E) {//Update Topological relation due to the reconnection of edge E
	//////////////////////For flipping of edge////////////////////////////////////
	int v1,v2,u1,u2;
	int e11,e12,e21,e22;
	int f1,f2,f3,f4;
	int Trial_Edge;
	bool State_Trial_Edge;//True: In Trial, False: conformed or rejected

	///////////////////////////////////0. Locate////////////////////////////////////
	//////substep1(in picture)
	v1=Vertices_of_Edge[E][0];
	v2=Vertices_of_Edge[E][1];

	//////substep2(in picture)
	u1=ExVertices_of_AdjFaces_of_Edge[E][0];
	u2=ExVertices_of_AdjFaces_of_Edge[E][1];
	f1=AdjFaces_of_Edge[E][0];
	f2=AdjFaces_of_Edge[E][1];
	int l;
	for(l=0;l<3;l++)if(Vertices_of_Face[f1][l]==v1)break;
	if(Vertices_of_Face[f1][(l+1)%3]!=v2){//if f1 is not v1->v2->xx, then define exchange f1 and f2
		swap(u1,u2);
		swap(f1,f2);
	}
	if(FindEdge(u1,u2)!=-1)return false;//{cout<<"Congratulations"<<endl;return false;}
	////pre-reject, in the project we move this segment to the other function, to make the performance better
	//if(Connection[v1].size()==4)return false;
	//if(Connection[v2].size()==4)return false;
	//if(Connection[u1].size()==8)return false;
	//if(Connection[u2].size()==8)return false;
	//////substep3(in picture)
	for(l=0;Vertices_of_Face[f1][l]!=v2;l++);
	e12=Edges_of_Face[f1][l];
	e11=Edges_of_Face[f1][(l+1)%3];

	for(l=0;Vertices_of_Face[f2][l]!=v1;l++);
	e21=Edges_of_Face[f2][l];
	e22=Edges_of_Face[f2][(l+1)%3];


	//////substep4(in picture)
	if(AdjFaces_of_Edge[e12][0]!=f1)f3=AdjFaces_of_Edge[e12][0];else f3=AdjFaces_of_Edge[e12][1];
	if(AdjFaces_of_Edge[e21][0]!=f2)f4=AdjFaces_of_Edge[e21][0];else f4=AdjFaces_of_Edge[e21][1];

	//////For debug
	//cout<<"Flipping:"<<E<<endl;
	//cout<<"u1="<<u1<<" u2="<<u2<<" v1="<<v1<<" v2="<<v2<<endl;
	//cout<<"f1="<<f1<<" f2="<<f2<<endl;
	//cout<<"e11="<<e11<<" e12="<<e12<<" e21="<<e21<<" e22="<<e22<<endl;
	///////////////////////////////////1. About E//////////////////////////////////
	Vertices_of_Edge[E][0]=u1;
	Vertices_of_Edge[E][1]=u2;
	if(AdjFaces_of_Edge[E][0]==f1) {
		ExVertices_of_AdjFaces_of_Edge[E][0]=v1;
		ExVertices_of_AdjFaces_of_Edge[E][1]=v2;
	}else {
		ExVertices_of_AdjFaces_of_Edge[E][0]=v2;
		ExVertices_of_AdjFaces_of_Edge[E][1]=v1;		
	}

	///////////////////////////////////2. About f1,f2/////////////////////////////
	for(int k=0;k<3;k++) {
		if(Vertices_of_Face[f1][k]==v2)Vertices_of_Face[f1][k]=u2;//make (v1,u1,u2)--->face f1
		if(Vertices_of_Face[f2][k]==v1)Vertices_of_Face[f2][k]=u1;//make (v2,u1,u2)--->face f2
	
		if(Edges_of_Face[f1][k]==E)Edges_of_Face[f1][k]=e21;
		if(Edges_of_Face[f1][k]==e12)Edges_of_Face[f1][k]=E;
		
		if(Edges_of_Face[f2][k]==E)Edges_of_Face[f2][k]=e12;
		if(Edges_of_Face[f2][k]==e21)Edges_of_Face[f2][k]=E;

		if(AdjFaces_of_Face[f1][k]==f2)AdjFaces_of_Face[f1][k]=f4;
		if(AdjFaces_of_Face[f1][k]==f3)AdjFaces_of_Face[f1][k]=f2;

		if(AdjFaces_of_Face[f2][k]==f1)AdjFaces_of_Face[f2][k]=f3;
		if(AdjFaces_of_Face[f2][k]==f4)AdjFaces_of_Face[f2][k]=f1;

	}
	


	///////////////////////////////////3. About u1,u2,v1,v2/////////////////////////
	//Modify connection
	auto it=Connection[v1].begin();
	for(;it<Connection[v1].end();it++)if((*it).y==E)break;
	Connection[v1].erase(it);

	for(it=Connection[v2].begin();it<Connection[v2].end();it++)if((*it).y==E)break;
	Connection[v2].erase(it);

	int2 tempInt2;
	tempInt2.x=u2;tempInt2.y=E;
	for(it=Connection[u1].begin();it<Connection[u1].end();it++)if((*it).y==e12)break;
	Connection[u1].insert(it,tempInt2);

	tempInt2.x=u1;tempInt2.y=E;
	for(it=Connection[u2].begin();it<Connection[u2].end();it++)if((*it).y==e21)break;
	Connection[u2].insert(it,tempInt2);

	//Modify Faces_of_Vertex
	auto It=Faces_of_Vertex[v1].begin();
	for(;It<Faces_of_Vertex[v1].end();It++)if((*It)==f2)break;
	Faces_of_Vertex[v1].erase(It);

	for(It=Faces_of_Vertex[v2].begin();It<Faces_of_Vertex[v2].end();It++)if((*It)==f1)break;
	Faces_of_Vertex[v2].erase(It);

	for(It=Faces_of_Vertex[u1].begin();It<Faces_of_Vertex[u1].end();It++)if((*It)==f3)break;
	Faces_of_Vertex[u1].insert(It,f2);

	for(It=Faces_of_Vertex[u2].begin();It<Faces_of_Vertex[u2].end();It++)if((*It)==f4)break;
	Faces_of_Vertex[u2].insert(It,f1);


	///////////////////////////////////4. About e11,e12,e21,e22/////////////////////
	if(AdjFaces_of_Edge[e21][0]==f2) {
		AdjFaces_of_Edge[e21][0]=f1;
		ExVertices_of_AdjFaces_of_Edge[e21][0]=u1;
	}else {
		AdjFaces_of_Edge[e21][1]=f1;
		ExVertices_of_AdjFaces_of_Edge[e21][1]=u1;
	}
	if(AdjFaces_of_Edge[e12][0]==f1) {
		AdjFaces_of_Edge[e12][0]=f2;
		ExVertices_of_AdjFaces_of_Edge[e12][0]=u2;
	}else {
		AdjFaces_of_Edge[e12][1]=f2;
		ExVertices_of_AdjFaces_of_Edge[e12][1]=u2;
	}

	if(ExVertices_of_AdjFaces_of_Edge[e11][0]==v2)ExVertices_of_AdjFaces_of_Edge[e11][0]=u2;else ExVertices_of_AdjFaces_of_Edge[e11][1]=u2;
	if(ExVertices_of_AdjFaces_of_Edge[e22][0]==v1)ExVertices_of_AdjFaces_of_Edge[e22][0]=u1;else ExVertices_of_AdjFaces_of_Edge[e22][1]=u1;


	///////////////////////////////////5. About f3,f4///////////////////////////////
	for(int k=0;k<3;k++) {
		if(AdjFaces_of_Face[f3][k]==f1)AdjFaces_of_Face[f3][k]=f2;
		if(AdjFaces_of_Face[f4][k]==f2)AdjFaces_of_Face[f4][k]=f1;
	}
	return true;
}

/*
void TriMesh::Print(int filenum) {
	cout<<"\n\n\n Print("<<filenum<<")"<<endl;
	cout<<"Number of vertices:"<<V<<" Number of Edges:"<<E<<" Number of Faces:"<<F<<endl;
	
	cout<<"Edges:"<<endl;
	if((AdjFaces_of_Edge.size()!=E)||(ExVertices_of_AdjFaces_of_Edge.size()!=E)||(Vertices_of_Edge.size()!=E))
		{cout<<"ERROR!"<<endl;return;}

	string EF("Edg");
	string TXT(".txt");
	string FF("Tri");
	string FNUM;
	FNUM=to_string(filenum);
	EF+=FNUM;
	EF+=TXT;
	ofstream Eout(EF);
	for(int k=0;k<E;k++) {
		cout<<"Edge "<<k<<" Vertices:";
		for(int l=0;l<Vertices_of_Edge[k].size();l++)cout<<Vertices_of_Edge[k][l]<<' ';

		cout<<" AdjFaces:";
		for(int l=0;l<AdjFaces_of_Edge[k].size();l++)cout<<AdjFaces_of_Edge[k][l]<<' ';
		cout<<" ExVerticles of Corresponding AdjFaces:";
		for(int l=0;l<ExVertices_of_AdjFaces_of_Edge[k].size();l++)cout<<ExVertices_of_AdjFaces_of_Edge[k][l]<<' ';
		cout<<endl;

		Eout<<k;
		for(int l=0;l<Vertices_of_Edge[k].size();l++)Eout<<' '<<Vertices_of_Edge[k][l];
		Eout<<endl;
	}
	Eout.close();

	FF+=FNUM;
	FF+=TXT;
	ofstream Fout(FF);

	cout<<"Faces:"<<endl;
	if((Vertices_of_Face.size()!=F)||(Edges_of_Face.size()!=F)||(AdjFaces_of_Face.size()!=F))
		{cout<<"ERROR!"<<endl;return;}

	for(int k=0;k<F;k++) {
		cout<<"Face "<<k<<" Vertices:";
		for(int l=0;l<Vertices_of_Face[k].size();l++)cout<<Vertices_of_Face[k][l]<<' ';
		Fout<<k;
		for(int l=0;l<Vertices_of_Face[k].size();l++)Fout<<' '<<Vertices_of_Face[k][l];
		Fout<<endl;
	
		cout<<" Edges:";
		for(int l=0;l<Edges_of_Face[k].size();l++)cout<<Edges_of_Face[k][l]<<' ';

		cout<<" AdjFaces:";
		for(int l=0;l<AdjFaces_of_Face[k].size();l++)cout<<AdjFaces_of_Face[k][l]<<' ';
		cout<<endl;
	}


	cout<<"Vertices:"<<endl;
	if((Faces_of_Vertex.size()!=V)||(Connection.size()!=V))
		{cout<<"ERROR!"<<endl;return;}

	for(int k=0;k<V;k++) {
		cout<<"Vertex "<<k<<" Faces:";
		for(int l=0;l<Faces_of_Vertex[k].size();l++)cout<<Faces_of_Vertex[k][l]<<' ';
	
		cout<<" Connections:(with, edge name):  ";
		for(int l=0;l<Connection[k].size();l++)cout<<"("<<Connection[k][l].x<<","<<Connection[k][l].y<<')';
		cout<<endl;
	}
}*/

bool TriMesh::Check() {
//////////////Step1: Check Faces/////////////
bool cor=true;
{
	int v1,v2,v3;
	int e1,e2,e3;
	int f1,f2,f3;
	for(int f=0;f<F;f++) {
		v1=Vertices_of_Face[f][0];v2=Vertices_of_Face[f][1];v3=Vertices_of_Face[f][2];
		e1=Edges_of_Face[f][0];e2=Edges_of_Face[f][1];e3=Edges_of_Face[f][2];
		f1=AdjFaces_of_Face[f][0];f2=AdjFaces_of_Face[f][1];f3=AdjFaces_of_Face[f][2];
		
		if(min(Vertices_of_Edge[e1][0],Vertices_of_Edge[e1][1])!=min(v1,v2)){cout<<"Face Error 0 at face"<<f<<endl;cor=false;}
		if(max(Vertices_of_Edge[e1][0],Vertices_of_Edge[e1][1])!=max(v1,v2)){cout<<"Face Error 1 at face"<<f<<endl;cor=false;}
		if(min(Vertices_of_Edge[e2][0],Vertices_of_Edge[e2][1])!=min(v2,v3)){cout<<"Face Error 2 at face"<<f<<endl;cor=false;}
		if(max(Vertices_of_Edge[e2][0],Vertices_of_Edge[e2][1])!=max(v2,v3)){cout<<"Face Error 3 at face"<<f<<endl;cor=false;}
		if(min(Vertices_of_Edge[e3][0],Vertices_of_Edge[e3][1])!=min(v3,v1)){cout<<"Face Error 4 at face"<<f<<endl;cor=false;}
		if(max(Vertices_of_Edge[e3][0],Vertices_of_Edge[e3][1])!=max(v3,v1)){cout<<"Face Error 5 at face"<<f<<endl;cor=false;}

		if((AdjFaces_of_Edge[e1][0]!=f)&&(AdjFaces_of_Edge[e1][1]!=f)){cout<<"Face Error 6 at face"<<f<<endl;cor=false;}
		if((AdjFaces_of_Edge[e1][0]!=f1)&&(AdjFaces_of_Edge[e1][1]!=f1)){cout<<"Face Error 7 at face"<<f<<endl;cor=false;}
		if((AdjFaces_of_Edge[e2][0]!=f)&&(AdjFaces_of_Edge[e2][1]!=f)){cout<<"Face Error 8 at face"<<f<<endl;cor=false;}
		if((AdjFaces_of_Edge[e2][0]!=f2)&&(AdjFaces_of_Edge[e2][1]!=f2)){cout<<"Face Error 9 at face"<<f<<endl;cor=false;}
		if((AdjFaces_of_Edge[e3][0]!=f)&&(AdjFaces_of_Edge[e3][1]!=f)){cout<<"Face Error 10 at face"<<f<<endl;cor=false;}
		if((AdjFaces_of_Edge[e3][0]!=f3)&&(AdjFaces_of_Edge[e3][1]!=f3)){cout<<"Face Error 11 at face"<<f<<endl;cor=false;}

		bool mark;
		mark=false;
		for(int k=0;k<Faces_of_Vertex[v1].size();k++)if(Faces_of_Vertex[v1][k]==f){mark=true;break;}
		if(!mark){cout<<"Face Error 12 at face"<<f<<endl;cor=false;}
		mark=false;
		for(int k=0;k<Faces_of_Vertex[v2].size();k++)if(Faces_of_Vertex[v2][k]==f){mark=true;break;}
		if(!mark){cout<<"Face Error 13 at face"<<f<<endl;cor=false;}
		mark=false;
		for(int k=0;k<Faces_of_Vertex[v3].size();k++)if(Faces_of_Vertex[v3][k]==f){mark=true;break;}
		if(!mark){cout<<"Face Error 14 at face"<<f<<endl;cor=false;}
	}
}
//////////////Step2: Check Edges/////////////
{
	int f1,f2;
	int v1,v2,u1,u2;
	bool mark;
	for(int e=0;e<E;e++) {
		f1=AdjFaces_of_Edge[e][0];
		f2=AdjFaces_of_Edge[e][1];
		v1=Vertices_of_Edge[e][0];
		v2=Vertices_of_Edge[e][1];
		u1=ExVertices_of_AdjFaces_of_Edge[e][0];
		u2=ExVertices_of_AdjFaces_of_Edge[e][1];
		mark=false;
		for(int l=0;l<=3;l++)if(Vertices_of_Face[f1][l]==u1){mark=true;break;}
		if(!mark){cout<<"Edge Error 0 at edge"<<e<<endl;cor=false;}
		mark=false;
		for(int l=0;l<=3;l++)if(Vertices_of_Face[f1][l]==v1){mark=true;break;}
		if(!mark){cout<<"Edge Error 1 at edge"<<e<<endl;cor=false;}
		mark=false;
		for(int l=0;l<=3;l++)if(Vertices_of_Face[f1][l]==v2){mark=true;break;}
		if(!mark){cout<<"Edge Error 2 at edge"<<e<<endl;cor=false;}

		mark=false;
		for(int l=0;l<=3;l++)if(Vertices_of_Face[f2][l]==u2){mark=true;break;}
		if(!mark){cout<<"Edge Error 3 at edge"<<e<<endl;cor=false;}
		mark=false;
		for(int l=0;l<=3;l++)if(Vertices_of_Face[f2][l]==v1){mark=true;break;}
		if(!mark){cout<<"Edge Error 4 at edge"<<e<<endl;cor=false;}
		mark=false;
		for(int l=0;l<=3;l++)if(Vertices_of_Face[f2][l]==v2){mark=true;break;}
		if(!mark){cout<<"Edge Error 5 at edge"<<e<<endl;cor=false;}

		int k1,k2;
		for(k1=0;k1<3;k1++)if(Vertices_of_Face[f1][k1]==u1)break;
		for(k2=0;k2<3;k2++)if(Vertices_of_Face[f2][k2]==u2)break;
		if((Vertices_of_Face[f1][(k1+1)%3]==v1)&&(Vertices_of_Face[f2][(k2+1)%3]!=v2)){cout<<"Edge Error 6 at edge"<<e<<endl;cor=false;}
		if((Vertices_of_Face[f1][(k1+1)%3]==v2)&&(Vertices_of_Face[f2][(k2+1)%3]!=v1)){cout<<"Edge Error 7 at edge"<<e<<endl;cor=false;}

		mark=false;
		for(int k=0;k<Connection[v1].size();k++)if((Connection[v1][k].x==v2)&&(Connection[v1][k].y==e)){mark=true;break;}
		if(!mark){cout<<"Edge Error 6 at edge"<<e<<endl;cor=false;}
		mark=false;
		for(int k=0;k<Connection[v2].size();k++)if((Connection[v2][k].x==v1)&&(Connection[v2][k].y==e)){mark=true;break;}
		if(!mark){cout<<"Edge Error 7 at edge"<<e<<endl;cor=false;}
	}
}
//////////////////step3 Check Vertices///////////////////////////////////////
{////检查连接关系中是否存在多余的信息
	int f;
	int2 c;
	int e;
	for(int v=0;v<V;v++) {
		for(int k=0;k<Faces_of_Vertex[v].size();k++) {
			f=Faces_of_Vertex[v][k];
			if((Vertices_of_Face[f][0]!=v)&&(Vertices_of_Face[f][1]!=v)&&(Vertices_of_Face[f][2]!=v)){cout<<"Vertex Error 0 at v"<<v<<endl;cor=false;}
			c=Connection[v][k];
			e=c.y;
			if(min(Vertices_of_Edge[e][0],Vertices_of_Edge[e][1])!=min(c.x,v)){cout<<"Vertex Error 1 at v"<<v<<endl;cor=false;}
			if(max(Vertices_of_Edge[e][0],Vertices_of_Edge[e][1])!=max(c.x,v)){cout<<"Vertex Error 2 at v"<<v<<endl;cor=false;}
		}

	}
}
{////检查时针顺序是否正确
	vector<int>e2s(0);
	for(int v=0;v<V;v++) {
		int f,f0;
		e2s.resize(0);
		f=Faces_of_Vertex[v][0];
		f0=f;
		int l;
		for(l=0;l<3;l++)if(Vertices_of_Face[f][l]==v)break;
		int v1,v2,e1,e2;
		v1=Vertices_of_Face[f][(l+1)%3];
		v2=Vertices_of_Face[f][(l+2)%3];
		e1=Edges_of_Face[f][l];
		e2=Edges_of_Face[f][(l+2)%3];
		do{
			e2s.push_back(e2);
			if(AdjFaces_of_Edge[e2][0]!=f)f=AdjFaces_of_Edge[e2][0];else f=AdjFaces_of_Edge[e2][1];
			for(l=0;l<3;l++)if(Vertices_of_Face[f][l]==v)break;
			if(Vertices_of_Face[f][(l+1)%3]!=v2){cout<<"Vertex Error 3 at v"<<v<<endl;cor=false;}
			v2=Vertices_of_Face[f][(l+2)%3];
			e2=Edges_of_Face[f][(l+2)%3];
		}while(f!=f0);
		for(l=0;l<Connection[v].size();l++)if(Connection[v][l].y==e2s[0])break;
		for(int k=0;k<Connection[v].size();k++)if(Connection[v][(l+k)%Connection[v].size()].y!=e2s[k]){cout<<"Vertex Error 4 at v"<<v<<endl;cor=false;}
	}
}

return cor;
}