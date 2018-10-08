#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<ctime>
#include <iomanip>
#include<cmath>

using namespace std;

#ifndef CELL_H
#define CELL_H

class Cell
{
	public:
		Cell(int no, float x, float y, float z);
		~Cell();
		int cellNo, collNo;
		float cellX ,cellY ,cellZ, velX, velY, velZ, porosity, regionVolume; 
		void increaseCollNo();
	private:

};
#endif
