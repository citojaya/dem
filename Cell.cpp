#include "Mill.h"


void Mill::setNeighbourCell(int cI, int nL[], int size)
{
	for (int i=0; i<size ; i++ )
	{
		cellCellList[cI][i] = nL[i];
	}
	cellCellListSize[cI] = size;
};


void Mill::deleteParticle(int cI, int pN)
{
	for (int i=0; i<cellParticleListSize[cI]; i++ )
	{
		if (cellParticleList[cI][i] == pN)
		{
			cellParticleList[cI][i] = cellParticleList[cI][cellParticleListSize[cI]-1];
			cellParticleListSize[cI] = cellParticleListSize[cI] - 1;
			break;
		}
	}
};

bool Mill::insertParticleToCell(int cI, int pN)
{
	cellParticleList[cI][cellParticleListSize[cI]] = pN;
	cellParticleListSize[cI] = cellParticleListSize[cI] + 1;
	if (cellParticleListSize[cI] > neighParListSize)
	{
		cout<<"particleListSize > "<<neighParListSize<<" (NeighbourCell)"<<endl;
		return false;
	}
	return true;
};
