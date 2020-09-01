/*
闭合曲面的三角网格中，三角形的以及点的时针方向需要保证。即便是在拓扑发生变化的情况下也是这样。
*/
#include "public.hpp"
class TriMesh {
	protected:

		///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////		
		///////////////////////////////////BASIC///////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////		
		int F,E,V;//Number of faces, number of edges, number of vertices




		///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////Topological characters of face//////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////
		vector<vector<int> >Vertices_of_Face;
		//Should conserve the clockwise
		
		vector<vector<int> >Edges_of_Face;
		/*
			Should conserve the clockwise
			Edges_of_Face[*][0] should be the edge with vertices(Vertices_of_Face[*][0],Vertices_of_Face[*][1])
			Edges_of_Face[*][1] should be the edge with vertices(Vertices_of_Face[*][1],Vertices_of_Face[*][2])
			Edges_of_Face[*][2] should be the edge with vertices(Vertices_of_Face[*][2],Vertices_of_Face[*][0])
		*/

		vector<vector<int> >AdjFaces_of_Face;
		/*
			Should conserve the clockwise
			The order of adjecent face is the same as their corresponding edge
		*/



		///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////Topological characters of edge////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////
		vector<vector<int> >Vertices_of_Edge;//no order

		vector<vector<int> >AdjFaces_of_Edge;
		vector<vector<int> >ExVertices_of_AdjFaces_of_Edge;
		//AdjFaces_of_Edge and ExVerticles_of_AdjFaces_of_Edge should have the same order
		
		//Topological characters of vertex
		vector<vector<int2> >Connection;//Should conserve the clockwise. *In an int2, the field x is the target of this connection, the field y is the name of the edge
		
		vector<vector<int> >Faces_of_Vertex;//Should conserve the clockwise
		//These two arrays don't need to keep the relative-order relation



		int FindOrCreate_Edge(int i,int j);
		int FindEdge(int i,int j);
	public:
		//In TriMesh.cpp
		TriMesh();
		//void Print(int);//For debug, Print all the topological relations
		bool Check();//For debug, check the consistency of topological relations
		bool Flip(int E);//Update Topological relation due to the reconnection of edge E
		int getV(){return this->V;}
		int getE(){return this->E;}
		int getF(){return this->F;}
};
