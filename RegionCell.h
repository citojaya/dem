#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<ctime>
#include <iomanip>
#include<cmath>
#include "Cell.h"
#define parSize 60
#define refCellSize 35

using namespace std;

#ifndef REGIONCELL_H
#define REGIONCELL_H

class RegionCell
{
	public:
		RegionCell(float x, float y, float z, float vol);
		~RegionCell();
		int collNo, noOfReferenceCells, noOfParticles;
		float cellX ,cellY ,cellZ; 
		float porosity, regionVolume;
		int particleList[parSize];
		Cell *refCellList[refCellSize];
		void insertRefCell(Cell *cl);
		void increaseCollNo();
	private:

};
#endif
