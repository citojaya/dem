/**
 * Neighbour cells 
 * author: Chandana Jayasundara 
 *			SIMPAS UNSW
 * ver: 4.0
 * 28 May 2005
 */
#include<cmath>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<ctime>
#define neighCellListSize 27
#define neighParListSize 20
//#include<omp.h>

using namespace std;

#ifndef NEIGHBOURCELL_H
#define NEIGHBOURCELL_H

class NeighbourCell
{
    public:
        NeighbourCell();
        ~NeighbourCell();
		int cellList[neighCellListSize]; // array for holding neighbour list cell numbers
		int particleList[neighParListSize]; //array for holding particle list within a cell
        void setNeighbourCell(int nL[], int size); 
		//void setCellParticleList(int pL[], int size); 
		void deleteParticle(int pN);
		bool insertParticle(int pN);
		int particleListSize,cellListSize;

    private:
		
        

};
#endif
