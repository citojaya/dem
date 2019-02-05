/*****************************************
 * Particle insertion and allocation 
 * author: Chandana Jayasundara 
 * SIMPAS UNSW
 * ver: 5.0
 * May 2007
 *****************************************/

#include "Mill.h"

/**************************** 
 * Program starts here
 ****************************/
void Mill::startMill()
{	
	shaftAngVel[2] = -2.0*M_PI*id->rotSpeed/60.0;
	shaftAngVel[0] = 0.0;
	shaftAngVel[1] = 0.0;
	rotorSpeed = -2.0*M_PI*id->rotSpeed/60.0;
	stepsBetweenFrames = (int)(2.0*M_PI/(noOfCfdFrames*fabs(rotorSpeed)*id->timeStep));
	
	//------ Read and update inintial CFD details ----------------
	// ifstream cfdFile("cfddetails.dat",ios::in);
	// int cfdCount = 0;
	// double cfxVX = 0.0;
	// double cfxVY = 0.0;
	// double cfxVZ = 0.0;
	// double cfxPX = 0.0;
	// double cfxPY = 0.0;
	// double cfxPZ = 0.0;
    //     if (!cfdFile)
    //     {
    //             cout<<"CFD INPUT FILE MISSING!!!!!!"<<endl;
    //             exit(1);
    //     }

	// while(cfdFile)
	// {
	// 	cfdFile.read(reinterpret_cast <char *>(&cfxVX),sizeof(double));
	// 	cfdFile.read(reinterpret_cast <char *>(&cfxVY),sizeof(double));
	// 	cfdFile.read(reinterpret_cast <char *>(&cfxVZ),sizeof(double));
	// 	cfdFile.read(reinterpret_cast <char *>(&cfxPX),sizeof(double));
	// 	cfdFile.read(reinterpret_cast <char *>(&cfxPY),sizeof(double));
	// 	cfdFile.read(reinterpret_cast <char *>(&cfxPZ),sizeof(double));
	// 	cfdCellVelX[cfdCount] = cfxVX;
	// 	cfdCellVelY[cfdCount] = cfxVY;
	// 	cfdCellVelZ[cfdCount] = cfxVZ;
	// 	cfdCellPGradX[cfdCount] = cfxPX;
	// 	cfdCellPGradY[cfdCount] = cfxPY;
	// 	cfdCellPGradZ[cfdCount] = cfxPZ;
	// 	cfdCount++;
	// }
	// cfdFile.close();
	//-----------------------------------------------------------
	// ifstream drumFile("in_drum",ios::in);
	// if (!drumFile)
	// {
	// 	cout<<"DRUM FILE MISSING!!!!!!"<<endl;
	// 	exit(1);
	// }
	// fstream File("out_shell.dat",ios::out | ios::binary);
	// fileCounter = 1;
	// ch = 's';
	
	// while(drumFile)
	// {
	// 	drumFile>>tempX>>tempY>>tempZ;
		
	// 	File.write(reinterpret_cast <const char *>(&tempX),sizeof(double));
	// 	File.write(reinterpret_cast <const char *>(&tempY),sizeof(double));
	// 	File.write(reinterpret_cast <const char *>(&tempZ),sizeof(double));
	// 	File.write(reinterpret_cast <const char *>(&fileCounter),sizeof(int));
	// }
		
	// drumFile.close();
	// fileCounter = 0;
	// ifstream shaftFile("in_shaft",ios::in);
	// if (!shaftFile)
	// {
	// 	cout<<"SHAFT FILE MISSING!!!!!!!"<<endl;
	// 	exit(1);
	// }

	// while(shaftFile)
	// {
	// 	if (fileCounter < 2644)
	// 	{
	// 		//finite element data points for the shaft//
	// 		shaftFile>>tempX>>tempY>>tempZ;
			
	// 		File.write(reinterpret_cast <const char *>(&tempX),sizeof(double));
	// 		File.write(reinterpret_cast <const char *>(&tempY),sizeof(double));
	// 		File.write(reinterpret_cast <const char *>(&tempZ),sizeof(double));
	// 		File.write(reinterpret_cast <const char *>(&fileCounter),sizeof(int));
	// 	}
	// 	else
	// 	{
	// 		shaftFile>>tempX>>tempY>>tempZ;
			
	// 		File.write(reinterpret_cast <const char *>(&tempX),sizeof(double));
	// 		File.write(reinterpret_cast <const char *>(&tempY),sizeof(double));
	// 		File.write(reinterpret_cast <const char *>(&tempZ),sizeof(double));
	// 	}
	// 	fileCounter++;
	// }
	// shaftFile.close();
	// File.close();

	ifstream infile("sample.dump",ios::in);
		
	//if sample.dump file does not exists//
	if(!infile)
	{
		infile.close();
		//ofstream dumpFile("sample.dump",ios::out);
                
        //while(index < 200)
	//		insertParticle();
		//insertParticle(0,19,5);
		//insertParticle(0,15.0,5.0);
	    	//insertParticle(1,-50,5);
	    	//insertParticle(0,-110,25);
                
		maxTime = 0.0;
		
	/*	while(index <id->noOfParticles)
		{
			maxTime = maxTime + 0.02*id->timeFactor;
            		counter = index + 5;
			while (index < counter)
			{
				if (index == id->noOfParticles)
				{
					break;
				}
				insertParticle();
				
			}
			cout<<" Particle "<<index<<" inserted"<<endl;
			allocateNeighbourList();
			//simulate();
		}
	*/
	}
	
	else
	{
		firstInsertion = false;
		maxTime = id->maxTime;
		char test = 'a';
		infile.ignore(1000,'#');
		infile>>test>>startingTime;
		infile.ignore(1000,'#');
		infile.ignore(1000,'#');
		while(infile)
		{
			infile>>test;
			if (test == '#')
			{
				break;
			}
			infile>>counter>>tempX>>tempY>>tempZ>>tempColor>>tempR>>dxDot>>dyDot>>dzDot>>xAngVel>>yAngVel>>zAngVel;
         
            //if (tempZ < 100)
            //{
			tempX = tempX*1e-3*id->lengthFactor;
			tempY = tempY*1e-3*id->lengthFactor;
			tempZ = tempZ*1e-3*id->lengthFactor;
			
			tempR = tempR*1e-3*id->lengthFactor;
			dxDot = dxDot*id->velocityFactor;
			dyDot = dyDot*id->velocityFactor;
			dzDot = dzDot*id->velocityFactor;

			xAngVel = xAngVel*id->angVelFactor;
			yAngVel = yAngVel*id->angVelFactor;
			zAngVel = zAngVel*id->angVelFactor;

			double tempRRR = sqrt(pow(tempX,2)+pow(tempY,2));
            if (allocateIntoCells(tempX,tempY,tempZ,noOfParticles))
            {
				particleNo[noOfParticles] = counter;
				//particleNo[noOfParticles] = noOfParticles;
				particleArrayIndex[noOfParticles] = noOfParticles;
				particleDiameter[noOfParticles] = tempR;
				particleDensity[noOfParticles] = id->pDensity;
				particleColor[noOfParticles] = tempColor;
				//particleMass[noOfParticles] = id->pMass;
                particleMass[noOfParticles] = (4.0/3.0)*M_PI*pow((0.5*tempR),3.0)*id->pDensity;
                particleInertia[noOfParticles] = 2*particleMass[noOfParticles]*pow(0.5*tempR,2)/5.0;

				particleX[noOfParticles] = tempX;
				particleY[noOfParticles] = tempY;
				particleZ[noOfParticles] = tempZ;

				particleVelX[noOfParticles] = dxDot;
				particleVelY[noOfParticles] = dyDot; 
				particleVelZ[noOfParticles] = dzDot;

				particleAngVelX[noOfParticles] = xAngVel;
				particleAngVelY[noOfParticles] = yAngVel;
				particleAngVelZ[noOfParticles] = zAngVel;
				noOfParticles++;
				index++;
            }
            //}
		}
		
		infile.ignore(1000,'#');
		infile>>test>>cfdFrameCounter;
		infile.close();
	}//end of else
/*
//------ Mirror image of sample.dump --------------------------
		ifstream infile2("sample2.dump",ios::in);
		infile2.ignore(1000,'#');
		infile2>>test>>startingTime;
		infile2.ignore(1000,'#');
		infile2.ignore(1000,'#');
		while(infile2)
		{
			infile2>>test;
			if (test == '#')
			{
				break;
			}
			infile2>>counter>>tempX>>tempY>>tempZ>>tempColor>>tempR>>dxDot>>dyDot>>dzDot>>xAngVel>>yAngVel>>zAngVel;
         
            //if (tempY < -5)
            //if (sqrt(tempX*tempX+tempY*tempY+tempZ*tempZ)<113.5)
	    //{
			tempX = tempX*1e-3*id->lengthFactor;
			tempY = tempY*1e-3*id->lengthFactor;
			tempZ = ((-(tempZ-102))+102)*1e-3*id->lengthFactor;
			
			tempR = tempR*1e-3*id->lengthFactor;
			dxDot = 0.0;//dxDot*id->velocityFactor;
			dyDot = 0.0;//dyDot*id->velocityFactor;
			dzDot = 0.0;//dzDot*id->velocityFactor;

			xAngVel = 0.0;//xAngVel*id->angVelFactor;
			yAngVel = 0.0;//yAngVel*id->angVelFactor;
			zAngVel = 0.0;//zAngVel*id->angVelFactor;

			double tempRRR = sqrt(pow(tempX,2)+pow(tempY,2));
            if (allocateIntoCells(tempX,tempY,tempZ,noOfParticles))
            {
				//particleNo[noOfParticles] = counter;
				particleNo[noOfParticles] = noOfParticles;
				particleArrayIndex[noOfParticles] = noOfParticles;
				particleDiameter[noOfParticles] = tempR;
				particleDensity[noOfParticles] = id->pDensity;
				particleColor[noOfParticles] = tempColor;
				//particleMass[noOfParticles] = id->pMass;
                                particleMass[noOfParticles] = (4.0/3.0)*M_PI*pow((0.5*tempR),3.0)*id->pDensity;
                                particleInertia[noOfParticles] = 2*particleMass[noOfParticles]*pow(0.5*tempR,2)/5.0;

				particleX[noOfParticles] = tempX;
				particleY[noOfParticles] = tempY;
				particleZ[noOfParticles] = tempZ;

				particleVelX[noOfParticles] = dxDot;
				particleVelY[noOfParticles] = dyDot; 
				particleVelZ[noOfParticles] = dzDot;

				particleAngVelX[noOfParticles] = xAngVel;
				particleAngVelY[noOfParticles] = yAngVel;
				particleAngVelZ[noOfParticles] = zAngVel;
				noOfParticles++;
				index++;
            }
            //}
		}
		
		infile2.ignore(1000,'#');
		infile2>>test>>rotorAngPosition;
		infile2.close();
*/
//-------------------------------------------------------------	
		//------ Set the angular positions of holes ------------
		for (int i=0; i<compIndex; i++)
		{
			if (componentType[i] == 1) //type is disk
			{
				for (int j=0; j<noOfSubComponents[i] ;j++ )
				{
					hollowCompTheta[i][j] = hollowCompTheta[i][j] - 2.0*M_PI*cfdFrameCounter/noOfCfdFrames;
				}
			}
		}

		//------ Set initial positions of disks -----------
		for(int i=0; i<diskFileArraySize1; i++)
		{
			//finite element points of the disk//
			tempX = pointArrayX[i];
			tempY = pointArrayY[i];
			tempR = sqrt(pow(tempX,2)+pow(tempY,2));
			if(tempX == 0 )
			{
				if (tempY > 0)
				{
					tempTheta = M_PI*0.5;
				}
				else
				{
					tempTheta = -M_PI*0.5;
				}

			}
			else if (tempX < 0)
			{
				tempTheta = atan(tempY/tempX) + M_PI;
			}
			else
			{
				tempTheta = atan(tempY/tempX);
			}
			tempTheta = tempTheta - M_PI/id->noOfHoles - 2.0*M_PI*cfdFrameCounter/noOfCfdFrames;
			/*
			if (tempTheta <= -2.0*M_PI)
			{
				tempTheta = tempTheta + 2.0*M_PI;
			}*/
			pointArrayX[i] = tempR*cos(tempTheta);
			pointArrayY[i] = tempR*sin(tempTheta);
		}
		//-------- End of disk position ------------------------

		cout<<" Particle "<<index<<" inserted"<<endl;
		//dumpFile.open("sample.dump",ios::out);
		maxTime = 0.0;
       
		/*while(index <id->noOfParticles)
		{
			counter = index + 2;
			while (index < counter)
			{
				if (index == id->noOfParticles)
				{
					break;
				}
				insertParticle();
			}
			cout<<" Particle "<<index<<" inserted"<<endl;
			cout<<" allocateNeighbourList()"<<endl;
			maxTime = maxTime + 0.03*id->timeFactor;
			allocateNeighbourList();
			//simulate();
		}*/
		cout<<"--- DONE  ---"<<endl;
		infile.close();
	//}//end of else
	allocateNeighbourList();
	maxTime = id->maxTime;
};

/*********************************************************
 * Insert particles into the system
 **********************************************************/
bool Mill::insertParticle()
//bool Mill::insertParticle(double x,double y,double z)
{
	dist =  id->shaftDia/2.0+id->particleDia/2.0;
	//generate a random number in between shaft diameter+particleDia/2 and 
	//chamberInnerDia-particleDia/2
	tempR = (rand()%(int)((id->chamberInnerDia-id->shaftDia)*0.5*1e3/id->lengthFactor -id->particleDia*1e3/id->lengthFactor-largestParDia)
										+((int)(id->shaftDia+id->particleDia)*0.5*1e3/id->lengthFactor + 1));
	//generate a random number in between 0 and 360
	tempTheta = (rand()%360+0);
	//generate a random number in between 0+particleDia/2  and chamberLength-particleDia/2
	tempZ = (rand()%(int)((id->chamberLength - id->particleDia)*1e3/id->lengthFactor - largestParDia)+(int)(id->particleDia*1e3*0.5/id->lengthFactor + 1));
	
	tempX = floor(tempR*cos(tempTheta*M_PI/180.0));
	tempZ = floor(tempZ);
	tempY = floor(tempR*sin(tempTheta*M_PI/180.0));

	tempX = tempX*1e-3*id->lengthFactor;
	tempY = tempY*1e-3*id->lengthFactor;
	tempZ = tempZ*1e-3*id->lengthFactor;
	//-- For single particle simulation ----
	//tempX = x*1e-3*id->lengthFactor;
	//tempY =	y*1e-3*id->lengthFactor;
	//tempZ =	z*1e-3*id->lengthFactor;
	//---------------------------------------

	//if((tempX !=pX || tempY!=pY || tempZ!=pZ) && tempY < -id->chamberInnerDia*0.33)
	//if((tempX !=pX || tempY!=pY || tempZ!=pZ) && tempY < -id->chamberInnerDia*0.25)
	//if((tempX !=pX || tempY!=pY || tempZ!=pZ) && tempY < -id->chamberInnerDia*0.125)
	//if((tempX !=pX || tempY!=pY || tempZ!=pZ) && tempY < id->chamberInnerDia*0)
	//if((tempX !=pX || tempY!=pY || tempZ!=pZ) && tempY < id->chamberInnerDia*0.125)
	if((tempX !=pX || tempY!=pY || tempZ!=pZ) && tempY < id->chamberInnerDia*0.5)
	//if((tempX !=pX || tempY!=pY || tempZ!=pZ) && tempY < 5)
	{
	//if(tempY > 0)
	//{
		pX = tempX;
		pY = tempY;
		pZ = tempZ;
		//set index=1 for the first time. Otherwise 'for' loop does not work for the first time
		if(firstInsertion)
		{
			xCoordinate = dist;
			zCoordinate = id->particleDia/2.0;
			yCoordinate = dist;
			if ((sqrt(pow(xCoordinate-pX,2)+
			pow(yCoordinate-pY,2)+
			pow(zCoordinate-pZ,2)) >= id->particleDia) &&
			sqrt(pow(pX,2)+pow(pY,2)) >= dist)
			{
				if(sqrt(pow(pX,2)+pow(pY,2))>=(id->diskDia+id->particleDia)/2.0 + 1*id->lengthFactor)
				{
					insertable = true;
				}
				else
				{
				for(int i=0; i<compIndex; i++)
				{
					if (componentType[i] == 1)//disk
					{
						//if((componentZ[i]-(id->particleDia+id->diskThick)/2.0 + 5*id->lengthFactor)>=pZ &&
						if((componentZ[i]-(id->particleDia+id->diskThick)/2.0)>=pZ &&
							pZ>=(componentZ[i]-id->diskApart+(-id->diskThick+id->particleDia)/2.0))
						{
							insertable = true;
							break;
						}
						//else if ((componentZ[i]+(id->particleDia+id->diskThick)/2.0 - 5*id->lengthFactor)<=pZ &&
						else if ((componentZ[i]+(id->particleDia+id->diskThick)/2.0)<=pZ &&
						pZ<=(componentZ[i]+id->diskApart+(id->diskThick-id->particleDia)/2.0))
						{
							insertable = true;
							break;
						}
						else
						{
							insertable = false;
							//break;
						}
					}
				}
				}//end of else
			}
			else
			{	
				insertable = false;
			}
		
		}
		else//if not first insertion
		{
			//current = totalParList->head;
			for (int i=0; i<noOfParticles; i++ )
			{
				xCoordinate = particleX[i];
				yCoordinate = particleY[i];
				zCoordinate = particleZ[i];

				if ((sqrt(pow(xCoordinate-pX,2)+
				pow(yCoordinate-pY,2)+
				pow(zCoordinate-pZ,2)) > (particleDiameter[i]+id->particleDia)/2.0) &&
				sqrt(pow(pX,2)+pow(pY,2)) > dist)
				{
					if(sqrt(pow(pX,2)+pow(pY,2))>(id->diskDia+id->particleDia)/2.0)
					{
						insertable = true;
					}
					else
					{
						for(int i=0; i<compIndex ; i++)
						{
							if (componentType[i] == 1)//disk
							{
								//if((componentZ[i]-(id->particleDia+id->diskThick)/2.0 + 5*id->lengthFactor)>=pZ &&
								if((componentZ[i]-(id->particleDia+id->diskThick)/2.0)>pZ &&
									pZ >(componentZ[i]-id->diskApart+(-id->diskThick+id->particleDia)/2.0))
								{
									insertable = true;
									break;
								}
								//else if ((componentZ[i]+(id->particleDia+id->diskThick)/2.0 - 5*id->lengthFactor)<=pZ &&
								else if ((componentZ[i]+(id->particleDia+id->diskThick)/2.0)<pZ &&
									pZ <(componentZ[i]+id->diskApart+(id->diskThick-id->particleDia)/2.0))
								{
									insertable = true;
									break;
								}
								else
								{
									insertable = false;
									//break;
								}
							}
						}
					}
				}
				else
				{
					insertable = false;
					break;
				}
			}//end of while
		}
		
		if(insertable)
		{
			if(firstInsertion)
			{
                if (allocateIntoCells(pX,pY,pZ,noOfParticles))
				{
					particleNo[noOfParticles] = index;
					particleX[noOfParticles] = pX;
					particleY[noOfParticles] = pY;
					particleZ[noOfParticles] = pZ;
			
					particleArrayIndex[noOfParticles] = noOfParticles;
					particleDiameter[noOfParticles] = id->particleDia;
					particleDensity[noOfParticles] = id->pDensity;
					particleMass[noOfParticles] = id->pMass;
					particleInertia[noOfParticles] = id->pInertia;
					particleColor[noOfParticles] = 1;
					firstInsertion = false;
					noOfParticles++;
					index++;
					return true;
				}
			}
			else
			{
				if (allocateIntoCells(pX,pY,pZ,noOfParticles))
				{
					particleNo[noOfParticles] = index;
					particleX[noOfParticles] = pX;
					particleY[noOfParticles] = pY;
					particleZ[noOfParticles] = pZ;
					particleArrayIndex[noOfParticles] = noOfParticles;
					particleDiameter[noOfParticles] = id->particleDia;
					particleDensity[noOfParticles] = id->pDensity;
					particleMass[noOfParticles] = id->pMass;
					particleInertia[noOfParticles] = id->pInertia;
					particleColor[noOfParticles] = 1;
					noOfParticles++;
					index++;
					return true;
				}
			}
			return true;
		}
		else
		return false;
	//}
	}
	return false;
};

/*******************************************************
 * At the begining each partilce is allocated to cells
 *******************************************************/
bool Mill::allocateIntoCells(double pX, double pY, double pZ, int pAIndex)
{
	for(int i=0; i<iIndex*jIndex*kIndex; i++)
	{
		if(fabs(cellCellX[i] - pX) <= cellDx/2.0)
        {
			if(fabs(cellCellY[i] - pY) <= cellDy/2.0)
			{
				if(fabs(cellCellZ[i] - pZ) <= cellDz/2.0)
				{	
					particleCellNo[pAIndex] = cellCellNo[i];
					if (!insertParticleToCell(i, pAIndex))
					{
                    }
					
                    particleInsertable[pAIndex] = true;
                    return true;
                }
			}
		}

	}//end of i loop
	return false;
};

/*********************************************
 * Allocate particles to the neighbour list
 *********************************************/
void Mill::allocateNeighbourList()
{
	for (int i=0; i<noOfParticles; i++ )
	{
        neighListSize[i] = 0; //reset neighbour list 
		int cellIndex = particleCellNo[i] - 1;
		
		if (cellIndex + 1 <= iIndex*jIndex) //All the cells where kIndex=0 (first Z layer)
		{
			for (int j=0; j<cellCellListSize[cellIndex]; j++ )//Neighbouring cells (max 27)
			{
				int neighCellIndex =  cellCellList[cellIndex][j] - 1; //array index of the neighbour cell
				for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
				{   
					int nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
					double xx = particleX[i] - particleX[nParIndex];
					double yy = particleY[i] - particleY[nParIndex];
					double zz = particleZ[i] - particleZ[nParIndex];
					if (neighCellIndex + 1 > iIndex*jIndex*(kIndex-1))// (last Z layer)
					{
						zz = particleZ[i] - (particleZ[nParIndex] - (id->diskThick + id->diskApart));
					}
						
					double parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
					if (parDistance <= rOut*id->lengthFactor)
					{
						insertNeighParticle(i,nParIndex);
					}
				}
			}
		}
		else if (cellIndex + 1 > iIndex*jIndex*(kIndex-1))//All the cells where kIndex=kIndex-1 (last Z layer)
		{
			for (int j=0; j<cellCellListSize[cellIndex]; j++ )//Neighbouring cells (max 27)
			{
				int neighCellIndex =  cellCellList[cellIndex][j] - 1; //array index of the neighbour cell
				for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
				{   
					int nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
					double xx = particleX[i] - particleX[nParIndex];
					double yy = particleY[i] - particleY[nParIndex];
					double zz = particleZ[i] - particleZ[nParIndex];
					if (neighCellIndex + 1 <= iIndex*jIndex)// (first Z layer)
					{
						zz = particleZ[i] - (particleZ[nParIndex]+ (id->diskThick + id->diskApart));
					}
						
					double parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
					if (parDistance <= rOut*id->lengthFactor)
					{
						insertNeighParticle(i,nParIndex);
					}
				}
			}
		}
		else 
		{
			for (int j=0; j<cellCellListSize[cellIndex]; j++ )//Neighbouring cells (max 27)
			{
				int neighCellIndex =  cellCellList[cellIndex][j] - 1; //array index of the neighbour cell
				for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
				{   
					int nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
					double xx = particleX[i] - particleX[nParIndex];
					double yy = particleY[i] - particleY[nParIndex];
					double zz = particleZ[i] - particleZ[nParIndex];
					double parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
					if (parDistance <= rOut*id->lengthFactor)
					{
						insertNeighParticle(i,nParIndex);
					}
				}
			}
		}
        //cout<<i<<" CELL INDEX  "<<cellIndex<<endl;
	}
};

/***************************************************
 * Each time step neighbour cells are updated
 ***************************************************/
void Mill::updateNeighList()
{
	int nParIndex = 0;
    int cellIndex = 0;
    int neighCellIndex = 0;
	double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xyz = 0.0;
    double parDistance = 0.0;
    int cellNum = 0;

	for (int i=0; i<noOfParticles; i++)
	{
		//if (particleY[i] > upperLimit)//***
		//{

		if (particleDisplacement[i] > allowedDistance*id->lengthFactor) //if loop 1
		{
		    particleDisplacement[i] = 0.0;
			cellNum = particleCellNo[i]; //cell number 
			if(fabs(cellCellX[cellNum-1] - particleX[i]) > cellDx/2.0 ||
						fabs(cellCellY[cellNum-1] - particleY[i]) > cellDy/2.0 ||
						fabs(cellCellZ[cellNum-1] - particleZ[i]) > cellDz/2.0) //if loop 6
				{
					deleteParticle((cellNum-1),particleArrayIndex[i]);
					
					for(int j=0; j<cellCellListSize[cellNum-1]; j++)
					{	
						int nCellIndex = cellCellList[cellNum-1][j]-1;
						if(particleInsertable[i])
						{
							if(fabs(cellCellX[nCellIndex] - particleX[i]) <= cellDx/2.0)
							{
								if(fabs(cellCellY[nCellIndex] - particleY[i])<= cellDy/2.0)
								{
									if(fabs(cellCellZ[nCellIndex] - particleZ[i]) <= cellDz/2.0)
									{
										particleCellNo[i] = cellCellNo[nCellIndex];
										particleInsertable[i] = false;
										insertParticleToCell(nCellIndex,particleArrayIndex[i]);
										break;
									}
								}
							}
								
						}//if particle insertable//

					}//end of j loop
				}//end of if loop 6

			//----------------------------------------------------------/

			//--------- Update the particles in the neighbourList -------
			for (int j=0; j<neighListSize[i]; j++ )
			{
				nParIndex = neighbourList[i][j];
				
				if (particleCellNo[i] <= iIndex*jIndex)//All the cells where kIndex=0 (first Z layer)
				{
					
					if (particleCellNo[nParIndex] > iIndex*jIndex*(kIndex-1))// (last Z layer)
					{
						xx = particleX[i] - particleX[nParIndex];
						yy = particleY[i] - particleY[nParIndex];
						zz = particleZ[i] - (particleZ[nParIndex] - (id->diskThick + id->diskApart));
						xyz = sqrt((xx*xx)+(yy*yy)+(zz*zz));
						if (xyz > rOut*id->lengthFactor)
						{
							deleteNeighParticle(nParIndex, i);
						}
					}
					else 
					{
						xx = particleX[i] - particleX[nParIndex];
						yy = particleY[i] - particleY[nParIndex];
						zz = particleZ[i] - particleZ[nParIndex];
						xyz = sqrt((xx*xx)+(yy*yy)+(zz*zz));
						if (xyz > rOut*id->lengthFactor)
						{
							deleteNeighParticle(nParIndex, i);
						}
					}
				}

				else if (particleCellNo[i] > iIndex*jIndex*(kIndex-1))//All the cells where kIndex=kIndex-1 (last Z layer)
				{
					nParIndex = neighbourList[i][j];
					if (particleCellNo[nParIndex] <= iIndex*jIndex)// (first Z layer)
					{
						xx = particleX[i] - particleX[nParIndex];
						yy = particleY[i] - particleY[nParIndex];
						zz = particleZ[i] - (particleZ[nParIndex] + (id->diskThick + id->diskApart));
						xyz = sqrt((xx*xx)+(yy*yy)+(zz*zz));
						if (xyz > rOut*id->lengthFactor)
						{
							deleteNeighParticle(nParIndex, i);
						}
					}
					else 
					{
						xx = particleX[i] - particleX[nParIndex];
						yy = particleY[i] - particleY[nParIndex];
						zz = particleZ[i] - particleZ[nParIndex];
						xyz = sqrt((xx*xx)+(yy*yy)+(zz*zz));
						if (xyz > rOut*id->lengthFactor)
						{
							deleteNeighParticle(nParIndex, i);
						}
					}
				}
				else 
				{
					nParIndex = neighbourList[i][j];
					xx = particleX[i] - particleX[nParIndex];
					yy = particleY[i] - particleY[nParIndex];
					zz = particleZ[i] - particleZ[nParIndex];
					xyz = sqrt((xx*xx)+(yy*yy)+(zz*zz));
					if (xyz > rOut*id->lengthFactor)
					{
						deleteNeighParticle(nParIndex, i);
					}
				}
			}
			//-----------------------------------------------------------------

			//------------ Update the new neighbourList -----------------------
			neighListSize[i] = 0; //set neighList size to zero in order to delete the neighbourList 
			cellIndex = particleCellNo[i] - 1;
			
			if (particleCellNo[i] <= iIndex*jIndex)//All the cells where kIndex=0 (first Z layer)
			{
				for (int j=0; j<cellCellListSize[cellIndex]; j++ )
				{
					neighCellIndex =  cellCellList[cellIndex][j] - 1; //array index of the neighbour cell
					if (neighCellIndex + 1 > iIndex*jIndex*(kIndex-1))// (last Z layer)
					{
						for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
						{
							nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
							xx = particleX[i] - particleX[nParIndex];
							yy = particleY[i] - particleY[nParIndex];
							zz = particleZ[i] - (particleZ[nParIndex] -(id->diskThick + id->diskApart));
							parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
							if (parDistance <= rOut*id->lengthFactor)
							{
								insertNeighParticle(i, nParIndex);
							}
						}
					}
					else 
					{
						for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
						{
							nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
							xx = particleX[i] - particleX[nParIndex];
							yy = particleY[i] - particleY[nParIndex];
							zz = particleZ[i] - particleZ[nParIndex];
							parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
							if (parDistance <= rOut*id->lengthFactor)
							{
								insertNeighParticle(i, nParIndex);
							}
						}

					}
				}
			}
			else if (particleCellNo[i] > iIndex*jIndex*(kIndex-1))//All the cells where kIndex=kIndex-1 (last Z layer)
			{
				for (int j=0; j<cellCellListSize[cellIndex]; j++ )
				{
					neighCellIndex =  cellCellList[cellIndex][j] - 1; //array index of the neighbour cell
					if (neighCellIndex + 1 <= iIndex*jIndex)// (first Z layer)
					{
						for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
						{
							nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
							xx = particleX[i] - particleX[nParIndex];
							yy = particleY[i] - particleY[nParIndex];
							zz = particleZ[i] - (particleZ[nParIndex] + (id->diskThick + id->diskApart));
							parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
							if (parDistance <= rOut*id->lengthFactor)
							{
								insertNeighParticle(i, nParIndex);
							}
						}
					}
					else 
					{
						for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
						{
							nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
							xx = particleX[i] - particleX[nParIndex];
							yy = particleY[i] - particleY[nParIndex];
							zz = particleZ[i] - particleZ[nParIndex];
							parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
							if (parDistance <= rOut*id->lengthFactor)
							{
								insertNeighParticle(i, nParIndex);
							}
						}

					}
				}
			}
			else 
			{
				for (int j=0; j<cellCellListSize[cellIndex]; j++ )
				{
					neighCellIndex =  cellCellList[cellIndex][j] - 1; //array index of the neighbour cell
					for (int k=0; k<cellParticleListSize[neighCellIndex]; k++)//for neighbour particle list
					{
						nParIndex = cellParticleList[neighCellIndex][k]; //array index of the particle in the neighbour cell
						xx = particleX[i] - particleX[nParIndex];
						yy = particleY[i] - particleY[nParIndex];
						zz = particleZ[i] - particleZ[nParIndex];
						parDistance = sqrt((xx*xx)+(yy*yy)+(zz*zz));
						if (parDistance <= rOut*id->lengthFactor)
						{
							insertNeighParticle(i, nParIndex);
						}
					}
				}
			}
			//-------------------------------------------------------------------------
		}// end of if particleDisplacement if loop 1
	//}//***
	}//end of i loop
};

/********************************************************
 * Simulation is carried out until total time = max time
 ********************************************************/
void Mill::simulate()
{
	int outFileCounter = 100;
	timeCounter = 0;
	int forceFileCounter = 0;
	//sampleDumpCounter/id->timeFactor > 0.005
	double sampleDumpCounter = 1.0/id->timeFactor;
	// if(cfdFrameCounter != 0)
	// {
	// 	updateCfdCells(cfdFrameCounter);//initial CFD cell update
	// }

	while(totalTime < maxTime)
	{
		for (int i=0;i<compIndex ;i++ )
		{
			for (int j=0;j<noOfSubComponents[i] ;j++ )
			{
				hollowCompTheta[i][j] = hollowCompTheta[i][j] + rotorSpeed*id->timeStep;
			}
		}

		//-------------Update the disk position ------------------------
		for(int i=0; i<diskFileArraySize1; i++)
		{
			//finite element points of the disk//
			tempX = pointArrayX[i];
			tempY = pointArrayY[i];
			tempR = sqrt(pow(tempX,2)+pow(tempY,2));
			if(tempX == 0 )
			{
				if (tempY > 0)
				{
					tempTheta = M_PI*0.5;
				}
				else
				{
					tempTheta = -M_PI*0.5;
				}

				//tempTheta = atan(tempY/tempX) + M_PI;
			}
			else if (tempX < 0)
			{
				tempTheta = atan(tempY/tempX) + M_PI;
			}
			else
			{
				tempTheta = atan(tempY/tempX);
			}
			
			tempTheta = tempTheta + rotorSpeed*id->timeStep;

			pointArrayX[i] = tempR*cos(tempTheta);
			pointArrayY[i] = tempR*sin(tempTheta);
		}
		//--------------- End of disk update ---------------------------------

		updateNeighList();
		startCalculation();
		
		rotorAngPosition = rotorAngPosition + rotorSpeed*id->timeStep;
		//diskAngPosition = diskAngPosition + rotorSpeed*id->timeStep;
        	cfdFrameAngle = cfdFrameAngle + rotorSpeed*id->timeStep;
	   	if (fabs(cfdFrameAngle) > 2.0*M_PI/noOfCfdFrames)
        	{
			cfdFrameCounter++;
			firstStep = false;
			stepNumber++;
			if (cfdFrameCounter > noOfCfdFrames)
			{
				cfdFrameCounter = 1; 
			}
	        	double simTime = totalTime/id->timeFactor;
            		double totTime = (totalTime+startingTime)/id->timeFactor;

            		writeDumpFile(simTime,totTime,cfdFrameCounter);
            		//writeCfdDetails();
			if (fabs(rotorAngPosition) > 2*M_PI)
			{
				rotorAngPosition = 0.0;	
		        	noOfCycles++;
				stepNumber = 1; //reset
                		cout<<"CYCLE "<<noOfCycles<<"  COMPLTED"<<endl;
            		}
			if(noOfCycles == id->setCycles)
            		{
                 		exit(0);
            		}
            		//updateCfdCells(cfdFrameCounter);
            		cfdFrameAngle = 0;
        	}
       
		//-------------------------------------------------
		//	THIS SECTION  FOR PARTICLE FORCE OUT FILE
		//-------------------------------------------------
        /*if (forceFileCounter >500) //every 0.02 seconds
		{
			//forceFile.seekp(0); //bring file pointer to the begining
			for (int i=0; i<noOfParticles; i++ )
			{
                double pX = particleX[i];
                double pY = particleY[i];
                double pZ = particleZ[i];
                double aCF = particleAverageCompForce[i];
                double aSF = particleAverageShearForce[i];
                
				//particle = totalParticleList[i];
				//arrayValue = particle->getCenter();
	            aveForceFile.write(reinterpret_cast <const char *>(&pX),sizeof(double));
	            aveForceFile.write(reinterpret_cast <const char *>(&pY),sizeof(double));
	            aveForceFile.write(reinterpret_cast <const char *>(&pZ),sizeof(double));
	            aveForceFile.write(reinterpret_cast <const char *>(&aCF),sizeof(double));
	            aveForceFile.write(reinterpret_cast <const char *>(&aSF),sizeof(double));
			}
			forceFileCounter = 0;
		}*/
		//---------- SECTION ENDS ----------------------------
 
		//if(outFileCounter >2)
		if(outFileCounter >250)
		//if(fabs(tempAngPosition) > (2.0*M_PI/id->noOfHoles - 0.05*M_PI/180.0))
		{
                        //writeDiskPosition();
                        writeParticlePosition();

                        for(int i=0; i<iIndex*jIndex*kIndex; i++)
                        {
                                double porosity = 0.0;
                                double totParVolume = 0.0;
                                for (int j=0; j<cellCellListSize[i]; j++)
                                {
                                   int neighCellIndex = cellCellList[i][j] - 1;
                                for (int m=0; m<cellParticleListSize[neighCellIndex]; m++)
                                {
                                        int parIndex = cellParticleList[neighCellIndex][m]; //particle index
                                        tempR = sqrt(pow((cellCellX[i] - particleX[parIndex]),2)
                                                                                        + pow((cellCellY[i] - particleY[parIndex]),2)
                                                                                        + pow((cellCellZ[i] - particleZ[parIndex]),2));
                                        if(tempR <= definedRegion*id->lengthFactor)//definedRegion is the radius of the spherical region
                                        {
                                                totParVolume  = totParVolume + findParVolume(particleDiameter[parIndex],definedRegion*id->lengthFactor,tempR);
                                        }
                                }
                                }
                                porosity = 1 - totParVolume/((4.0/3.0)*pow(definedRegion*id->lengthFactor*1.0,3));
                cellPorosity[i] = porosity;
                        }

                        //writeZonePosition();
            outFileCounter = 0;
                        //tempAngPosition = 0;
			zoneCounter++;
		}

        moveParticle();
		totalTime = totalTime + id->timeStep;
		outFileCounter++;
		tempAngPosition = tempAngPosition + rotorSpeed*id->timeStep;
		forceFileCounter++;
			
		if(timeCounter > 49)
		{
			cout<<" TIME "<<totalTime/id->timeFactor<<endl;
			timeCounter = 0;
		}
		timeCounter++;
		
	}//end while//
	ch = 'e';
	//File.write(reinterpret_cast <const char *>(&ch),sizeof(char));

	//----- After maximum time calculate wear data near hole region ----------
	for (int i=0;i<noOfSubComponents[0] ;i++ ) //no of holes
	{
		double holeAngle = hollowCompTheta[0][i];
		double holeCenterVec[3]; //vector representing hole centre
		holeCenterVec[0] = hollowCompR[0][i]*cos(holeAngle);
		holeCenterVec[1] = hollowCompR[0][i]*sin(holeAngle);
		holeCenterVec[2] = 0.0;//hollowCompZ[0][i];
		double holeEdgeRadiusVec[3];//vector representing radius to hole edge
		holeEdgeRadiusVec[0] = (hollowCompR[0][i] + hollowCompDia[0][i]*0.5)*cos(holeAngle);
		holeEdgeRadiusVec[1] = (hollowCompR[0][i] + hollowCompDia[0][i]*0.5)*sin(holeAngle);
		holeEdgeRadiusVec[2] = 0.0;//hollowCompZ[0][i];
		double holeRadiusVec[3]; //vector representing hole radius
		holeRadiusVec[0] = holeEdgeRadiusVec[0] - holeCenterVec[0];
		holeRadiusVec[1] = holeEdgeRadiusVec[1] - holeCenterVec[1];
		holeRadiusVec[2] = 0.0;//holeEdgeRadiusVec[2] - holeCenterVec[2];
		double rotatedAngle = 2.0*M_PI/noOfHoleDivisions; 
		for (int j=0; j<noOfHoleDivisions; j++ )
		{
			double r_holeRadiusVec[3]; //rotated hole radius vector
			r_holeRadiusVec[0] = holeRadiusVec[0]*cos(rotatedAngle*j) + holeRadiusVec[1]*sin(rotatedAngle*j);
			r_holeRadiusVec[1] = holeRadiusVec[1]*cos(rotatedAngle*j) - holeRadiusVec[0]*sin(rotatedAngle*j);
			r_holeRadiusVec[2] = holeRadiusVec[2];
			double resultantVec[3]; //vector representing a point on hole edge
			resultantVec[0] = holeCenterVec[0] + r_holeRadiusVec[0];
			resultantVec[1] = holeCenterVec[1] + r_holeRadiusVec[1];
			resultantVec[2] = hollowCompZ[0][i] + 0.5*hollowCompLength[0][i];
			for (int k=0; k<diskFileArraySize1; k++)
			{
				if (sqrt(pow((resultantVec[0] - pointArrayX[k]),2)+pow((resultantVec[1] - pointArrayY[k]),2)+pow((resultantVec[2] - pointArrayZ[k]),2)) < wearGapLimit*id->lengthFactor)
				{
					holeWearData1[j] = (1.0*holeWearData1[j]*k)/(k+1) + (1.0/(k+1))*pointData1[k]; 
					holeWearData2[j] = (1.0*holeWearData2[j]*k)/(k+1) + (1.0/(k+1))*pointData2[k]; 
				}
			}//end of k loop
			errorFile<<setw(12)<<holeWearData1[j]<<setw(12)<<holeWearData2[j]<<endl;
			holeWearData1[j] = 0.0;
			holeWearData2[j] = 0.0;
		}//end of j loop

		errorFile<<endl;

	}//end of i loop
	//------------------------------------------------------------------------
	errorFile<<"totalWear1  "<<totalWear1<<endl;
	errorFile<<"totalWear2  "<<totalWear2<<endl;
};

double Mill::findParVolume(double pD, double defR, double d)
{
    double ro = pD*0.5;
    double Ro = defR;
    double volume = 0.0;
    if(d <= Ro - ro)
    {
		volume = (4.0/3.0)*pow(ro,3);
    }
    else if(d >= Ro+ro)
    {
		volume = 0.0;
    }
    else
    {
		double L2 = ro;
		double L1 = 0.5*(pow(ro,2)+pow(d,2)-pow(Ro,2))/d;
		double L3 = 0.5*(pow(Ro,2)+pow(d,2)-pow(ro,2))/d;
		double L4 = Ro;
		double H2 = L2*pow(ro,2)-pow(L2,3)/3.0;
		double H1 = L1*pow(ro,2)-pow(L1,3)/3.0;

		double H22 = L4*pow(Ro,2)-pow(L4,3)/3.0;
		double H11 = L3*pow(Ro,2)-pow(L3,3)/3.0;
		volume = H2-H1+H22-H11;
    }
	return volume;	
};
void Mill::writeCfdDetails()
{
	fstream cfdFile("cfddetails.dat",ios::out | ios::binary);
	for(int i=0; i<cfd_iIndex*cfd_jIndex*cfd_kIndex; i++)
	{
		double cfdVelX = cfdCellVelX[i];
		double cfdVelY = cfdCellVelY[i];
		double cfdVelZ = cfdCellVelZ[i];
		double cfdPX = cfdCellPGradX[i];
		double cfdPY = cfdCellPGradY[i];
		double cfdPZ = cfdCellPGradZ[i];
		cfdFile.write(reinterpret_cast <const char *>(&cfdVelX),sizeof(double));
		cfdFile.write(reinterpret_cast <const char *>(&cfdVelY),sizeof(double));
		cfdFile.write(reinterpret_cast <const char *>(&cfdVelZ),sizeof(double));
		cfdFile.write(reinterpret_cast <const char *>(&cfdPX),sizeof(double));
		cfdFile.write(reinterpret_cast <const char *>(&cfdPY),sizeof(double));
		cfdFile.write(reinterpret_cast <const char *>(&cfdPZ),sizeof(double));
	}
	cfdFile.close();
};
/*******************************************************************
 * Simulation details such as particle position, particle velocity,
 * disk postions are written to a text file at equal time intervals
 *******************************************************************/
 
void Mill::writeDumpFile(double sT, double tT,int cfdFC)
{
	ofstream dumpFile("sample.dump",ios::out);//set file pointer at the begining
	dumpFile<<"SIMULATION TIME "<<sT<<endl;
	//dumpFile<<"SIMULATION TIME "<<totalTime<<endl;
	dumpFile<<"TIME "<<endl;
	dumpFile<<"#"<<endl;
	//dumpFile<<"+ "<<totalTime+startingTime<<endl;
	dumpFile<<"+ "<<tT<<endl;

	dumpFile<<"#"<<endl;
	
	dumpFile<<"np "<<noOfParticles<<endl;
	dumpFile<<endl;
	dumpFile<<"   Particle Position  "<<endl; 
	dumpFile<<"#"<<endl;

	for (int i=0; i<noOfParticles; i++ )
	{
		double tempXX = particleX[i]*1e3/id->lengthFactor;
		double tempYY = particleY[i]*1e3/id->lengthFactor;
		double tempZZ = particleZ[i]*1e3/id->lengthFactor;
	
		double dxxDot = particleVelX[i]/id->velocityFactor;
		double dyyDot =	particleVelY[i]/id->velocityFactor;
		double dzzDot =	particleVelZ[i]/id->velocityFactor;
		double xxAngVel = particleAngVelX[i]/id->angVelFactor;  
		double yyAngVel = particleAngVelY[i]/id->angVelFactor;  
		double zzAngVel = particleAngVelZ[i]/id->angVelFactor;  
		dumpFile<<"+ "<<setw(8)<<fixed<<setprecision(4)<<particleNo[i]
				<<setw(10)<<tempXX<<setw(10)<<tempYY<<setw(10)<<tempZZ
				<<setw(3)<<particleColor[i]<<setw(10)<<particleDiameter[i]*1e3/id->lengthFactor
				<<setw(12)<<dxxDot<<setw(12)<<dyyDot<<setw(12)<<dzzDot
				<<setw(12)<<xxAngVel<<setw(12)<<yyAngVel<<setw(12)<<zzAngVel<<endl;	
	}
	
	dumpFile<<"#"<<endl;
	
	dumpFile<<"   Rotor angular position /(rad)"<<endl;
	dumpFile<<"#"<<endl;
	dumpFile<<"+ "<<cfdFC<<endl;
	//dumpFile<<"+ "<<rotorAngPosition/id->freqFactor<<endl;
	dumpFile<<"#"<<endl;
	dumpFile.close();
};


void Mill::writeDiskPosition()
{
	ch = 's';
    int rank = 0;
    char firstStr[] = "out_";
    char midStr[] = "disk";
    char lastStr[] = ".dat";
    char rankStr[3];
    sprintf(rankStr,"%d",rank);
    strcat(firstStr,rankStr);
    strcat(firstStr,midStr);
    char numStr[3];
    sprintf(numStr,"%d",zoneCounter);
    strcat(firstStr,numStr);
    strcat(firstStr,lastStr);
    fstream File(firstStr,ios::out | ios::binary);

	for(int i=0; i<diskFileArraySize1; i++)
	{
		//finite element points of the disk//
		double tempXX = pointArrayX[i]/id->lengthFactor;
		double tempYY = pointArrayY[i]/id->lengthFactor;

		File.write(reinterpret_cast <const char *>(&tempXX),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempYY),sizeof(double));
		int tempInt = 1+i;
		double tempVal = pointArrayZ[i]/id->lengthFactor;
		double pData1 = pointData1[i];
		double pData2 = pointData2[i];
		File.write(reinterpret_cast <const char *>(&tempVal),sizeof(double));
		File.write(reinterpret_cast <const char *>(&pData1),sizeof(double));
		File.write(reinterpret_cast <const char *>(&pData2),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempInt),sizeof(int));
	}

	for(int i=0; i<diskFileArraySize2; i++)
	{
		tempX = elementArrayX[i];
		tempY = elementArrayY[i];
		tempZ = elementArrayZ[i];
		File.write(reinterpret_cast <const char *>(&tempX),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempY),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempZ),sizeof(double));
	}
	File.close();
	
};


void Mill::writeParticlePosition()
{
    dataFile<<" TIME = "<<totalTime/id->timeFactor<<endl;
	for (int i=0; i<noOfParticles ;i++ )
	{
		dataFile<<fixed<<setprecision(3)<<setw(10)<<particleX[i]*1e3/id->lengthFactor<<setw(10)
		<<particleY[i]*1e3/id->lengthFactor
		<<setw(10)<<particleZ[i]*1e3/id->lengthFactor
		<<setw(10)<<particleVelX[i]/id->velocityFactor
		<<setw(10)<<particleVelY[i]/id->velocityFactor
		<<setw(10)<<particleVelZ[i]/id->velocityFactor
		<<setw(10)<<particleDiameter[i]*1e3/id->lengthFactor<<endl;
	}
	dataFile<<endl;
	
}

/*
void Mill::writeParticlePosition()
{
    int rank = 0;
    char firstStr[] = "out_";
    char midStr[] = "par";
    char lastStr[] = ".dat";
    char rankStr[3];
    sprintf(rankStr,"%d",rank);
    strcat(firstStr,rankStr);
    strcat(firstStr,midStr);
    char numStr[3];
    sprintf(numStr,"%d",zoneCounter);
    strcat(firstStr,numStr);
    strcat(firstStr,lastStr);
    fstream File(firstStr,ios::out | ios::binary);

	double totalKE = 0.0;//total kinetic energy of particles
	//totalPE = 0.0;//total potential energy of particles
	double tempVel = 0.0;
	//double tempShear = 0.0;
	//double tempKE = 0.0;
	//double tempCompF = 0.0;
	double tempDia = 0.0;
	int tempNoColl = 0;
	int parNo = 0;
	int parCol = 0;// color of the particle
	ch = 's';
	//outputFile<<" PARTICLE "<<endl;
	for (int i=0; i<noOfParticles; i++ )
	{
		tempX = particleX[i]/id->lengthFactor; //mm
		tempY = particleY[i]/id->lengthFactor; //mm
		tempZ = particleZ[i]/id->lengthFactor; //mm
		tempDia = particleDiameter[i]/id->lengthFactor;
		parNo = particleNo[i];
		parCol = particleColor[i];
		tempNoColl = totalCollisions[i];
		dxDot = particleVelX[i]/id->velocityFactor;//(m/s) 
		dyDot = particleVelY[i]/id->velocityFactor;//(m/s) 
		dzDot = particleVelZ[i]/id->velocityFactor;//(m/s) 
		tempVel = floor((1e+4)*(pow(dxDot,2)+pow(dyDot,2)+pow(dzDot,2)));
		tempVel = tempVel*(1e-4);
		totalKE = totalKE + tempVel;
		double dFX = particleDragForceX[i];
		double dFY = particleDragForceY[i];
		double dFZ = particleDragForceZ[i];
		double dPX = particlePForceX[i];
		double dPY = particlePForceY[i];
		double dPZ = particlePForceZ[i];
		File.write(reinterpret_cast <const char *>(&parNo),sizeof(int));
		File.write(reinterpret_cast <const char *>(&tempX),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempY),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempZ),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dxDot),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dyDot),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dzDot),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempDia),sizeof(double));
		File.write(reinterpret_cast <const char *>(&parCol),sizeof(int));
		File.write(reinterpret_cast <const char *>(&tempNoColl),sizeof(int));
		File.write(reinterpret_cast <const char *>(&dFX),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dFY),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dFZ),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dPX),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dPY),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dPZ),sizeof(double));
	}
	File.close();
};
*/
/*
void Mill::writePositionData()
{
	for (int i=0; i<noOfParticles; i++ )
	{
		int parNo = particleNo[i]; 
        double aCF = particleAverageCompForce[i]; 
		double aSF = particleAverageShearForce[i];
		double radialVel = 0.0;
		double tangentialVel = 0.0;
		
		tempX = particleX[i];
		tempY = particleY[i];
		tempZ = particleZ[i];
		double tempR = sqrt(pow(tempX,2)+pow(tempY,2));
		double theta = 0.0;
		if (tempX != 0)
		{
			theta = atan(fabs(tempY/tempX));
		}
		else if (tempX == 0)
		{
			theta = M_PI*0.5;
		}
		dxDot = particleVelX[i]*(1e-3);//(m/s) 
		dyDot = particleVelY[i]*(1e-3);//(m/s) 
		dzDot = particleVelZ[i]*(1e-3);//(m/s) 
		if (tempY >= 0)
		{
			if (tempX >= 0)
			{
				radialVel = dxDot*cos(theta) + dyDot*sin(theta);
				tangentialVel = dxDot*sin(theta) - dyDot*cos(theta);
			}
			else 
			{
				radialVel = -dxDot*cos(theta) + dyDot*sin(theta);
				tangentialVel = dxDot*sin(theta) + dyDot*cos(theta);
			}
		}
		
		else 
		{
			if (tempX >= 0)
			{
				radialVel = dxDot*cos(theta) - dyDot*sin(theta);
				tangentialVel = -dxDot*sin(theta) - dyDot*cos(theta);
			}
			else 
			{
				radialVel = -dxDot*cos(theta) - dyDot*sin(theta);
				tangentialVel = -dxDot*sin(theta) + dyDot*cos(theta);
			}
		}

		particleDataFile.write(reinterpret_cast <const char *>(&parNo),sizeof(int));
		particleDataFile.write(reinterpret_cast <const char *>(&tempX),sizeof(double));
		particleDataFile.write(reinterpret_cast <const char *>(&tempY),sizeof(double));
		particleDataFile.write(reinterpret_cast <const char *>(&tempZ),sizeof(double));
		particleDataFile.write(reinterpret_cast <const char *>(&tempR),sizeof(double));
		particleDataFile.write(reinterpret_cast <const char *>(&aCF),sizeof(double));
		particleDataFile.write(reinterpret_cast <const char *>(&aSF),sizeof(double));
		particleDataFile.write(reinterpret_cast <const char *>(&radialVel),sizeof(double));
		particleDataFile.write(reinterpret_cast <const char *>(&tangentialVel),sizeof(double));
	}
}; 
*/
void Mill::writeZonePosition()
{
	ch = 's';
    int rank = 0;
    char firstStr[] = "out_";
    char midStr[] = "zone";
    char lastStr[] = ".dat";
    char rankStr[3];
    sprintf(rankStr,"%d",rank);
    strcat(firstStr,rankStr);
    strcat(firstStr,midStr);
    char numStr[3];
    sprintf(numStr,"%d",zoneCounter);
    strcat(firstStr,numStr);
    strcat(firstStr,lastStr);
    fstream File(firstStr,ios::out | ios::binary);

	for(int i=0; i<iIndex*jIndex*kIndex; i++)
	{
		double totalForceX = 0.0;
		double totalForceY = 0.0;
		double totalForceZ = 0.0;
		double reynoldsNo = 0.0;
		dxDot = 0.0;
		dyDot = 0.0;
		dzDot = 0.0;
		double shearForce = 0.0;
		double vel = 0.0;
		double tempVal = 0.0;
        double collE = 0.0;
		int aveNoOfCollisions = 0;
		tempX = cellCellX[i]/id->lengthFactor; 
		tempY = cellCellY[i]/id->lengthFactor; 
		tempZ = cellCellZ[i]/id->lengthFactor; 
		double iE = cellImpactEnergy[i];
		double sE = cellShearEnergy[i];
		double kE = cellCollisionEnergy[i];
		double tkE = cellTangCollisionEnergy[i];
		double porosity = cellPorosity[i];
		for (int m=0; m<cellParticleListSize[i]; m++)
		{	
			int parArrayIndex = cellParticleList[i][m];
			//shearForce = shearForce + particleAverageShearForce[parArrayIndex]/id->forceFactor;
			//compForce = compForce + particleAverageCompForce[parArrayIndex]/id->forceFactor;
			//collE = collE + particleAverageCollE[parArrayIndex]/id->energyFactor;

			dxDot = dxDot + particleVelX[parArrayIndex]; 
			dyDot = dyDot + particleVelY[parArrayIndex]; 
			dzDot = dzDot + particleVelZ[parArrayIndex]; 
			totalForceX = totalForceX + particleDragForceX[parArrayIndex];// + particlePForceX[parArrayIndex];
			totalForceY = totalForceY + particleDragForceY[parArrayIndex];// + particlePForceY[parArrayIndex];
			totalForceZ = totalForceZ + particleDragForceZ[parArrayIndex];// + particlePForceZ[parArrayIndex];
			reynoldsNo = reynoldsNo + particleReynolds[parArrayIndex];
		}
	
				
		if (cellParticleListSize[i]!= 0)
		{
			dxDot = dxDot/(cellParticleListSize[i]*id->velocityFactor);
			dyDot = dyDot/(cellParticleListSize[i]*id->velocityFactor);
			dzDot = dzDot/(cellParticleListSize[i]*id->velocityFactor);
			totalForceX = totalForceX*id->volumeFactor/(id->forceFactor*cellVolume)/porosity;
			totalForceY = totalForceY*id->volumeFactor/(id->forceFactor*cellVolume)/porosity;
			totalForceZ = totalForceZ*id->volumeFactor/(id->forceFactor*cellVolume)/porosity;
			reynoldsNo = reynoldsNo/cellParticleListSize[i];
		}

		int noOfCollisions = cellNoOfContacts[i];
		int noOfPar = cellParticleListSize[i];
		File.write(reinterpret_cast <const char *>(&tempX),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempY),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tempZ),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dxDot),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dyDot),sizeof(double));
		File.write(reinterpret_cast <const char *>(&dzDot),sizeof(double));
		File.write(reinterpret_cast <const char *>(&totalForceX),sizeof(double));
		File.write(reinterpret_cast <const char *>(&totalForceY),sizeof(double));
		File.write(reinterpret_cast <const char *>(&totalForceZ),sizeof(double));
		File.write(reinterpret_cast <const char *>(&iE),sizeof(double));
		File.write(reinterpret_cast <const char *>(&sE),sizeof(double));
		File.write(reinterpret_cast <const char *>(&kE),sizeof(double));
		File.write(reinterpret_cast <const char *>(&tkE),sizeof(double));
		File.write(reinterpret_cast <const char *>(&reynoldsNo),sizeof(double));
		File.write(reinterpret_cast <const char *>(&porosity),sizeof(double));
		File.write(reinterpret_cast <const char *>(&noOfPar),sizeof(int));
		File.write(reinterpret_cast <const char *>(&noOfCollisions),sizeof(int));
	}
	File.close();
};


/*******************************************************
 * Set the next position and velocity of each particle
 *******************************************************/
 void Mill::moveParticle()
{
	for (int i=0; i<noOfParticles; i++ )
	{
		//if (particleY[i] > upperLimit)//***
		//{

                particleInsertable[i] = true;
                //particleAverageCompForce[i] = 0.0;
                //particleAverageShearForce[i] = 0.0;

                dxDot = (particleForceX[i]*id->timeStep)/particleMass[i];
                dyDot = (particleForceY[i]*id->timeStep)/particleMass[i];
                dzDot = (particleForceZ[i]*id->timeStep)/particleMass[i];

                particleVelX[i] = particleVelX[i] + dxDot;
                particleVelY[i] = particleVelY[i] + dyDot;
                particleVelZ[i] = particleVelZ[i] + dzDot;

                xAngVel = particleAngVelX[i] + (particleMomentumX[i])*id->timeStep/particleInertia[i];
                yAngVel = particleAngVelY[i] + (particleMomentumY[i])*id->timeStep/particleInertia[i];
                zAngVel = particleAngVelZ[i] + (particleMomentumZ[i])*id->timeStep/particleInertia[i];

                particleAngVelX[i] = xAngVel;
                particleAngVelY[i] = yAngVel;
                particleAngVelZ[i] = zAngVel;

		double dX = particleVelX[i]*id->timeStep;
		double dY = particleVelY[i]*id->timeStep;
		double dZ = particleVelZ[i]*id->timeStep;
		particleDisplacement[i] = particleDisplacement[i] + sqrt((dX*dX)+(dY*dY)+(dZ*dZ));
		particleX[i] = particleX[i] + dX;
		particleY[i] = particleY[i] + dY;
		particleZ[i] = particleZ[i] + dZ;
		
		//--- With boundary conditions ----------------------------------------
		/*if (particleZ[i] < 0)
		{
			particleZ[i] = particleZ[i] + id->chamberLength;
			particleDisplacement[i] = 2*allowedDistance*id->lengthFactor;
		}
		if (particleZ[i] > id->chamberLength)
		{
			particleZ[i] = particleZ[i] - id->chamberLength;
			particleDisplacement[i] = 2*allowedDistance*id->lengthFactor;
		}*/
		//---------------------------------------------------------------------
		
		particleForceX[i]= 0.0;    
		particleForceY[i]= 0.0;    
		particleForceZ[i]= 0.0;    
		particleMomentumX[i]= 0.0; 
		particleMomentumY[i]= 0.0; 
		particleMomentumZ[i]= 0.0; 
		deleteContactList(i);
		
		//}//***
	}//end of particle list//
};


void Mill::deleteContactList(int parI)
{
	instantContactListSize[parI] = 0;
};

void Mill::setContactList(int pI, int nParI)
{
    instantContactList[pI][instantContactListSize[pI]] = nParI;
    instantContactListSize[pI] = instantContactListSize[pI] + 1;
	if (instantContactListSize[pI] > collisionIndex)
	{
		cout<<" instantContactListSize[pI] > "<<collisionIndex<<endl;
		exit(0);
	}
};


void Mill::insertNeighParticle(int parI, int pI)
{
	neighbourList[parI][neighListSize[parI]] = pI;
	neighListSize[parI] = neighListSize[parI] + 1;
    if(neighListSize[parI] > neighListArraySz)
    {
        cout<<" neighListSize > "<<neighListArraySz<<endl;
        exit(0);
    }
};

void Mill::deleteNeighParticle(int parI, int pI)
{
	for (int i=0; i<neighListSize[parI] ;i++ )
	{
		if (neighbourList[parI][i] == pI)
		{
			neighbourList[parI][i] = neighbourList[parI][neighListSize[parI]-1];
			neighListSize[parI] = neighListSize[parI] - 1;
			break;
		}
	}
};



