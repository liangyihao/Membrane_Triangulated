//deal with hard repulsions
//range:[-L/2,L/2]
//hard wall boundary condition
#include "public.hpp"
#define MaxNumberParticlePerCell 15
typedef struct CellStruct{
	int num_particle;//number of particles in this cell
	int particle[MaxNumberParticlePerCell];//index of particles in this cell

	CellStruct(){
		num_particle=0;
	}
	void Insert(int i) {
		if(num_particle==MaxNumberParticlePerCell){cout<<"Insertion failed. Need a larger cell volume"<<endl;return;}
		particle[num_particle]=i;
		num_particle++;
		return;
	}
	void Delete(int i) {
		for(int k=0;k<num_particle;k++)
			if(particle[k]==i) {
				particle[k]=particle[num_particle-1];
				num_particle--;
				return;
			}
		cout<<"Deletion failed. The particle "<<i<<" not found."<<endl;return;
	}
}Cell;

class CellList{
protected:
	vector<double3>X;
	double L,radius,dL;//size of each cell, should be equals 2*radius
	int NC;//Number of cells per dimension
	vector<int3>InWhichCell;
	vector<Cell>Cells;//CellList(i*NC*NC+j*NC+k) access the cell[i][j][k]
public:
	CellList(double L,double radius);
	bool Overlap(int I,double3 NX);
	void Update(int I,double3 NX);
};


