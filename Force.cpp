#include "Mill.h"

/**************************************************************************
 * for each iteration this method is called by ParticleMotion
 * For the total particle list finds the contact particle and if gap between
 * particle exists calculate the contact force. Find the contact surface for 
 * each particle and if overlap exists find the contact force
 **************************************************************************/
void Mill::startCalculation()
{	
	torque = 0.0;
    dsmaxCff = id->sfricp*(2.0 - id->pPoisonR)/(2*(1.0 - id->pPoisonR));
	torqueOnMill[0] = 0.0;
	torqueOnMill[1] = 0.0;
	torqueOnMill[2] = 0.0;

	for (int i=0; i<noOfParticles; i++ )
	{
		//if (particleY[i] > upperLimit)//***
		//{
	    dragFX = 0.0;
        dragFY = 0.0;
        dragFZ = 0.0;

        parX = particleX[i];
        parY = particleY[i];
        parZ = particleZ[i];
        parDia = particleDiameter[i];
		num = particleCellNo[i];

		collEnergy = 0.0;
		instImpactEnergy = 0.0;
		parIndex = i; //current particle
        particleDragForceX[i] =  0.0;//dragFX;
        particleDragForceY[i] =  0.0;//dragFY;
        particleDragForceZ[i] =  0.0;//dragFZ;
		
		//------ Setting up existing contacts ---------------------------  
		
		for (int k=0; k<instantContactListSize[parIndex]; k++ )
		{
			int arrIndex = instantContactList[parIndex][k];//arrayIndex of the contact particle
			particleContactCalculated[arrIndex] = true;
		}
		//-----------------------------------------------------------------
        //------- CFD code --------------------------------------------------
	
		//------------Drag force due to the liquid --------------------------
        int velCellNum = num;
		double cellVX = 0.0;//cellArrayList[velCellNum-1]->velX*id->velocityFactor;
		double cellVY = 0.0;//cellArrayList[velCellNum-1]->velY*id->velocityFactor;
		double cellVZ = 0.0;//cellArrayList[velCellNum-1]->velZ*id->velocityFactor;
		double cellPressureGradX = 0.0;
		double cellPressureGradY = 0.0;
		double cellPressureGradZ = 0.0;

		//double parNewX = parX*cos(-rotorAngPosition) - parY*sin(-rotorAngPosition);
		//double parNewY = parX*sin(-rotorAngPosition) + parY*cos(-rotorAngPosition);
        double parNewX = parX;
		double parNewY = parY;
		parNewX = parNewX + 0.5*id->chamberInnerDia;
			parNewY = parNewY + 0.5*id->chamberInnerDia;
			int ii = (int)floor(parNewX/cfdDx);
		if(ii < 0)
		{
		   ii = 0;
		}
		else if(ii > cfd_iIndex-1)
		{
		   ii = cfd_iIndex-1;
		}
		int jj = (int)floor(parNewY/cfdDy);
		if(jj < 0)
		{
		   jj = 0;
		}
		else if(jj > cfd_jIndex-1)
		{
			jj = cfd_jIndex-1;
        }
        int kk = (int)floor(parZ/cfdDz);
        if(kk < 0)
        {
           kk =0;
        }
        else if(kk > cfd_kIndex-1)
        {
            kk = cfd_kIndex-1;
        }
        int p = ii;
        int q = jj*cfd_jIndex;
        int r = kk*cfd_iIndex*cfd_jIndex;
        int tempCellIndex = p+q+r;
        if (tempCellIndex < 0 || tempCellIndex >cfd_iIndex*cfd_jIndex*cfd_kIndex)
        {
            cout<<" x  "<<parNewX/id->lengthFactor<<" y  "<<parNewY/id->lengthFactor<<"  z  "<<parZ/id->lengthFactor<<endl;
			cout<<p<<"  "<<q<<"  "<<r<<endl;
			cout<<" tempCellIndex "<<tempCellIndex<<endl;
			exit(0);
        }


		/*cellVX = (cellVelX[tempCellIndex]*cos(rotorAngPosition) - cellVelY[tempCellIndex]*sin(rotorAngPosition))*id->velocityFactor;  
		cellVY = (cellVelX[tempCellIndex]*sin(rotorAngPosition) + cellVelY[tempCellIndex]*cos(rotorAngPosition))*id->velocityFactor;  
		cellVZ = cellVelZ[tempCellIndex]*id->velocityFactor;  
		cellPressureGradX = (cellPGradX[tempCellIndex]*cos(rotorAngPosition) - cellPGradY[tempCellIndex]*sin(rotorAngPosition))*id->pressureFactor; 
		cellPressureGradY = (cellPGradX[tempCellIndex]*sin(rotorAngPosition) + cellPGradY[tempCellIndex]*cos(rotorAngPosition))*id->pressureFactor;
		cellPressureGradZ = cellPGradZ[tempCellIndex]*id->pressureFactor; 
		*/
		//-- Standard equations to find out drag coefficient -------------------
/*
		if(firstStep)
		{
			cellVX = cfdCellVelX[tempCellIndex]*id->velocityFactor;
			cellVY = cfdCellVelY[tempCellIndex]*id->velocityFactor;
			cellVZ = cfdCellVelZ[tempCellIndex]*id->velocityFactor;
	                cellPressureGradX = cfdCellPGradX[tempCellIndex]*id->pressureFactor;
                	cellPressureGradY = cfdCellPGradY[tempCellIndex]*id->pressureFactor;
                	cellPressureGradZ = cfdCellPGradZ[tempCellIndex]*id->pressureFactor;

		}
		else
		{*/
			cellVX =(pre_cfdCellVelX[tempCellIndex] 
					+  (cfdCellVelX[tempCellIndex] 
					- pre_cfdCellVelX[tempCellIndex])*frameStepCounter/stepsBetweenFrames)
					*id->velocityFactor;
                        cellVY =(pre_cfdCellVelY[tempCellIndex]
                                        +  (cfdCellVelY[tempCellIndex]
                                        - pre_cfdCellVelY[tempCellIndex])*frameStepCounter/stepsBetweenFrames)
                                        *id->velocityFactor;
                        cellVZ =(pre_cfdCellVelZ[tempCellIndex]
                                        +  (cfdCellVelZ[tempCellIndex]
                                        - pre_cfdCellVelZ[tempCellIndex])*frameStepCounter/stepsBetweenFrames)
                                        *id->velocityFactor;
	                cellPressureGradX = (pre_cfdCellPGradX[tempCellIndex]
                                        +  (cfdCellPGradX[tempCellIndex]
                                        - pre_cfdCellPGradX[tempCellIndex])*frameStepCounter/stepsBetweenFrames)
                                        *id->pressureFactor;
                        cellPressureGradY = (pre_cfdCellPGradY[tempCellIndex]
                                        +  (cfdCellPGradY[tempCellIndex]
                                        - pre_cfdCellPGradY[tempCellIndex])*frameStepCounter/stepsBetweenFrames)
                                        *id->pressureFactor;
                        cellPressureGradZ = (pre_cfdCellPGradZ[tempCellIndex]
                                        +  (cfdCellPGradZ[tempCellIndex]
                                        - pre_cfdCellPGradZ[tempCellIndex])*frameStepCounter/stepsBetweenFrames)
                                        *id->pressureFactor;

                //}

		double totParVolume = 0.0;
		for(int j=0; j<neighListSize[parIndex]; j++)
		{
			tempR = sqrt(pow((particleX[parIndex] - particleX[neighbourList[parIndex][j]]),2)
							+ pow((particleY[parIndex] - particleY[neighbourList[parIndex][j]]),2)
								+ pow((particleZ[parIndex] - particleZ[neighbourList[parIndex][j]]),2));
			if(tempR <= definedRegion*id->lengthFactor)//definedRegion is the radius of the spherical region
			{
				totParVolume  = totParVolume + findParVolume(particleDiameter[parIndex],definedRegion*id->lengthFactor,tempR);
			}
		}
		double instPorosity = 1.0 - totParVolume/((4.0/3.0)*pow(definedRegion*id->lengthFactor*1.0,3));

		double coeff = 0.0;
		double dragCoeff = 0.0;
		double modifiedPorosity = 0.0;
		double relativeVel = sqrt((cellVX-particleVelX[i])*(cellVX-particleVelX[i]) + (cellVY-particleVelY[i])*(cellVY-particleVelY[i]) + (cellVZ-particleVelZ[i])*(cellVZ-particleVelZ[i]));
		double reynoldsNo = (instPorosity*parDia/id->lengthFactor*(relativeVel/id->velocityFactor)*fluidDensity)/(fluidViscosity);
		//-- Standard equations -------------------------------------------------
		if (reynoldsNo < 1e-6)//if Reynolds number is too small following parameters cannot be defined
		{
			coeff = 0.0;
			dragCoeff = 0.0;
			modifiedPorosity = 0.0;
		}
		else 
		{
			coeff = 3.7-0.65*exp(-(1.5-log10(reynoldsNo))*(1.5-log10(reynoldsNo))*0.5);
			dragCoeff = (0.63+4.8/sqrt(reynoldsNo))*(0.63+4.8/sqrt(reynoldsNo));
			modifiedPorosity = pow(instPorosity,(-coeff));
		}
		//------------------------------------------------------------------------
		
		dragFX = 0.5*modifiedPorosity*dragCoeff*fluidDensity*id->densityFactor*instPorosity*(cellVX - particleVelX[i])*(fabs(cellVX - particleVelX[i]))*M_PI*pow(parDia*0.5,2)
					+ (-cellPressureGradX*M_PI*pow(particleDiameter[parIndex],3))/6.0;
        dragFY = 0.5*modifiedPorosity*dragCoeff*fluidDensity*id->densityFactor*instPorosity*(cellVY - particleVelY[i])*(fabs(cellVY - particleVelY[i]))*M_PI*pow(parDia*0.5,2)
					+ (-cellPressureGradY*M_PI*pow(particleDiameter[parIndex],3))/6.0;
        dragFZ = 0.5*modifiedPorosity*dragCoeff*fluidDensity*id->densityFactor*instPorosity*(cellVZ - particleVelZ[i])*(fabs(cellVZ - particleVelZ[i]))*M_PI*pow(parDia*0.5,2)
				+ (-cellPressureGradZ*M_PI*pow(particleDiameter[parIndex],3))/6.0;
		bouyancyF = M_PI*pow(parDia,3)*fluidDensity*id->densityFactor/6.0;
		particleDragForceX[i] =  0.5*modifiedPorosity*dragCoeff*fluidDensity*id->densityFactor*instPorosity*(cellVX - particleVelX[i])*(fabs(cellVX - particleVelX[i]))*M_PI*pow(parDia*0.5,2);
		particleDragForceY[i] =  0.5*modifiedPorosity*dragCoeff*fluidDensity*id->densityFactor*instPorosity*(cellVY - particleVelY[i])*(fabs(cellVY - particleVelY[i]))*M_PI*pow(parDia*0.5,2);
		particleDragForceZ[i] = 0.5*modifiedPorosity*dragCoeff*fluidDensity*id->densityFactor*instPorosity*(cellVZ - particleVelZ[i])*(fabs(cellVZ - particleVelZ[i]))*M_PI*pow(parDia*0.5,2);
		particlePForceX[i] = -cellPressureGradX*M_PI*pow(particleDiameter[parIndex],3)/6.0;
		particlePForceY[i] = -cellPressureGradY*M_PI*pow(particleDiameter[parIndex],3)/6.0;
		particlePForceZ[i] = -cellPressureGradZ*M_PI*pow(particleDiameter[parIndex],3)/6.0;
		//---- end of CFD code -----------------------------------------

		parX = particleX[i];     
		parY = particleY[i]; 
		parZ = particleZ[i]; 
		parDia = particleDiameter[i];
		num = particleCellNo[i];

		pForce[0] =  particleForceX[i];
		pForce[1] =  particleForceY[i];
		pForce[2] =  particleForceZ[i];
		pMom[0] =  particleMomentumX[i];    
		pMom[1] =  particleMomentumY[i];
		pMom[2] =  particleMomentumZ[i];

		//-- Initialise ------------
		momentum[0] = 0.0;
		momentum[1] = 0.0;
		momentum[2] = 0.0;
		momentumR[0] = 0.0;
		momentumR[1] = 0.0;
		momentumR[2] = 0.0;
		rotorForce[0] = 0.0;
		rotorForce[1] = 0.0;
		rotorForce[2] = 0.0;
		//---------------------------

		gForce = particleMass[i];		
		pForce[0] = pForce[0] + dragFX;
		pForce[1] = pForce[1] + dragFY + bouyancyF + (-gForce);
		//pForce[1] = pForce[1] + dragFY;
		pForce[2] = pForce[2] + dragFZ;
		
		particleForceX[i] = pForce[0];
		particleForceY[i] = pForce[1];
		particleForceZ[i] = pForce[2];
		if(num < 1 || num > iIndex*jIndex*kIndex)
		{
			cout<< "'ERROR' ARRAY INDEX OUT OF BOUNDS-Force.cpp line 162 num="<<num<<endl;;
			cout<<"PAR "<<particleNo[i]<<"  "<<parX<<" "<<parY<<" "<<parZ<<endl;
			exit(1);
		}
		/*
		//----- With boundary conditions -------------------------------------------------------------------
		if (num <= iIndex*jIndex)//All the cells where kIndex=0 (first Z layer)
		{
			for (int k=0; k<neighListSize[i]; k++)//for neighbour particle list
			{
				nParIndex = neighbourList[i][k]; //array index of the particle in the neighbour cell
				if(particleNo[nParIndex] != particleNo[i]) // To make sure 'nPar' is not 'par' //
				{
					if (particleContactCalculated[nParIndex] == false) // to check the contact calculation being done 
					{
						if (particleCellNo[nParIndex] > iIndex*jIndex*(kIndex-1))// (last Z layer) 
						{                                                                          
							nParZ = particleZ[nParIndex]- (id->diskThick + id->diskApart);
						}
						else
						{
							nParZ = particleZ[nParIndex];
						}

						gap = sqrt(((parX- particleX[nParIndex])*(parX- particleX[nParIndex]))+
										((parY- particleY[nParIndex])*(parY- particleY[nParIndex]))+
										((parZ- nParZ)*(parZ- nParZ)))-(parDia + particleDiameter[nParIndex])*0.5;
						collisionGap = gap;
						if (collisionGap > collisionGapCoeff)
						{
							deleteContact(parIndex, particleArrayIndex[nParIndex]);
							deleteContact(nParIndex, particleArrayIndex[parIndex]);
						}
						if(gap < 0)
						{ 
							calForceBetweenPart();
							setContactList(nParIndex,parIndex);
							setContactList(parIndex,nParIndex);
						}
					}
				}
			}
		}
		
		else if (num  > iIndex*jIndex*(kIndex-1))//All the cells where kIndex=kIndex-1 (last Z layer)
		{
			for (int k=0; k<neighListSize[i]; k++)//for neighbour particle list
			{
				nParIndex = neighbourList[i][k]; //array index of the particle in the neighbour cell
				if(particleNo[nParIndex] != particleNo[i])  // To make sure 'nPar' is not 'par' //
				{
					if (particleContactCalculated[nParIndex] == false)// to check the contact calculation being done 
					{
						if (particleCellNo[nParIndex] <= iIndex*jIndex)// (first Z layer) 
						{                                                                 
							nParZ = particleZ[nParIndex] + (id->diskThick + id->diskApart);
						}
						else
						{
							nParZ = particleZ[nParIndex];
						}
						gap = sqrt(((parX- particleX[nParIndex])*(parX- particleX[nParIndex]))+
										((parY- particleY[nParIndex])*(parY- particleY[nParIndex]))+
										((parZ- nParZ)*(parZ- nParZ)))-(parDia + particleDiameter[nParIndex])*0.5;
						collisionGap = gap;
						if (collisionGap > collisionGapCoeff)
						{
							deleteContact(parIndex, particleArrayIndex[nParIndex]);
							deleteContact(nParIndex, particleArrayIndex[parIndex]);
						}
						if(gap < 0)
						{ 
							calForceBetweenPart();
							setContactList(nParIndex,parIndex);
							setContactList(parIndex,nParIndex);
						}
					}
				}
			}
		}
		
		else
		{*/
		//----------------------------------------------------------------------------------------------------
		//------- Without boundary conditions -----------------------------------------------------------------
			for (int k=0; k<neighListSize[i]; k++)//for neighbour particle list
			{
				nParIndex = neighbourList[i][k]; //array index of the particle in the neighbour cell
				if(particleNo[nParIndex] != particleNo[i]) // To make sure 'nPar' is not 'par' //
				{
					if (particleContactCalculated[nParIndex] == false)// to check the contact calculation being done 
					{
						nParZ = particleZ[nParIndex];
						gap = sqrt(((parX- particleX[nParIndex])*(parX- particleX[nParIndex]))+
										((parY- particleY[nParIndex])*(parY- particleY[nParIndex]))+
										((parZ- nParZ)*(parZ- nParZ)))-(parDia + particleDiameter[nParIndex])*0.5;
						collisionGap = gap;
						if (collisionGap > collisionGapCoeff)
						{
							deleteContact(parIndex, particleArrayIndex[nParIndex]);
							deleteContact(nParIndex, particleArrayIndex[parIndex]);
						}
						if(gap < 0)
						{ 
							calForceBetweenPart();
							setContactList(nParIndex,parIndex);
							setContactList(parIndex,nParIndex);
						}
					}
				}

			}//end of k loop (neighbour particle list)
		//}
		//reset "particleContactCalculated" before going to the next particle
		//-------------------------------------------------------------------------------------------------------
		
		//reset "particleContactCalculated" before going to the next particle
		for (int k=0; k<instantContactListSize[parIndex]; k++ )
		{
			int arrIndex = instantContactList[parIndex][k];//arrayIndex of the contact particle
			particleContactCalculated[arrIndex] = false;
		}
        findContact();
		if (particleNoOfContacts[i] != 0)
		{
			//particleAverageCompForce[i] = particleAverageCompForce[i]/particleNoOfContacts[i]; 
			//particleAverageShearForce[i] = particleAverageShearForce[i]/particleNoOfContacts[i];
		    //particleAverageCollE[i] = particleAverageCollE[i]/particleNoOfContacts[i];
        }
		//File<<endl;
	}//end of while for particle list//
	//}//***
	totalTorque = totalTorque + torque;
	if(timeCount > 80)
	{	
		torqueFile<<totTime/id->timeFactor<<"  "<<totalTorque/(id->momentFactor*(timeCount+1))<<endl;
        timeCount = 0;
		totalTorque = 0.0;
	}

	totTime = totTime + id->timeStep;
	timeCount++;
	timeCount1++;
	frameStepCounter++;
	/*if (fabs(cfdFrameAngle) > 2.0*M_PI/noOfCfdFrames)
	{
		cout<<"MEMORY DELETED"<<endl;
		delete []cfxCellVelX;
		delete []cfxCellVelY;
		delete []cfxCellVelZ;
		delete []cfxCellPressure;
		delete []cfdCellNeighbourDistance;
        delete []cfdCellNeighbourPressure;
        delete []cfdCellNeighbourVelX;
        delete []cfdCellNeighbourVelY;
        delete []cfdCellNeighbourVelZ;
		delete []tempCfxCellX;
		delete []tempCfxCellY;
		cfdFrameAngle = 0;
	}*/
};

void Mill::updateCfdCells(int fcounter)
{
	cout<<"CFD frame counter "<<fcounter<<endl;
	cout<<"step counter  "<<stepNumber<<endl;
	cout<<"rotorAngPosition "<<rotorAngPosition<<endl;
	cout<<"noOfCycles "<<noOfCycles<<endl;
	/*for (int i=0; i<cfd_iIndex*cfd_jIndex*cfd_kIndex; i++ )
	{
		cfdCellVelX[i] = 0.0;
		cfdCellVelY[i] = 0.0;
		cfdCellVelZ[i] = 0.0;
		cfdCellPressure[i] = 0.0;
	}*/
	frameStepCounter = 0;
	int cfd_counter = 0;		 
	double  cfd_angleRotated =(2.0*M_PI/noOfCfdFrames)*(fcounter);
	double cfx_cellX = 0.0;
	double cfx_cellY = 0.0;
	double cfx_cellZ = 0.0;
	double cfx_cellVelX = 0.0;
	double cfx_cellVelY = 0.0;
	double cfx_cellVelZ = 0.0;
	double cfx_cellPressure = 0.0;
    char firstStr[] = "Step-";
    char fileNoStr[3];
    //fileNoStr[0] = ' ';
    //fileNoStr[1] = ' ';
    //fileNoStr[2] = ' ';
    sprintf(fileNoStr,"%d",fcounter);
    strcat(firstStr,fileNoStr);
	
	ifstream cfdFile(firstStr,ios::in);
	if (!cfdFile)
	{
		cout<<"Step-"<<fcounter<<"  missing !!!!!!!!!!!!!"<<endl;
		exit(0);
	}
	while (cfdFile)
	{
		cfdFile>>cfx_cellPressure>>cfx_cellVelX>>cfx_cellVelY>>cfx_cellVelZ;
		
		tempCfxCellX[cfd_counter] = cfxCellX[cfd_counter]*cos(cfd_angleRotated) + cfxCellY[cfd_counter]*sin(cfd_angleRotated);   
		tempCfxCellY[cfd_counter] = -cfxCellX[cfd_counter]*sin(cfd_angleRotated) + cfxCellY[cfd_counter]*cos(cfd_angleRotated);   
		cfxCellPressure[cfd_counter] = cfx_cellPressure;
		cfxCellVelX[cfd_counter] = cfx_cellVelX*cos(cfd_angleRotated) + cfx_cellVelY*sin(cfd_angleRotated);
		cfxCellVelY[cfd_counter] = -cfx_cellVelX*sin(cfd_angleRotated) + cfx_cellVelY*cos(cfd_angleRotated);
		
		cfxCellVelZ[cfd_counter] = cfx_cellVelZ;
		cfd_counter++;
	}

	for (int m=0; m<cfd_counter ;m++ )
	{
		//exit(0);
		//-- shift X and Y axes to eliminate minus (-) coordinates
		cfx_cellX = tempCfxCellX[m] + 0.5*id->chamberInnerDia;
		cfx_cellY = tempCfxCellY[m] + 0.5*id->chamberInnerDia;
		cfx_cellZ = cfxCellZ[m]+ 0.5*id->chamberLength;
		cfx_cellVelX = cfxCellVelX[m];
		cfx_cellVelY = cfxCellVelY[m];
		cfx_cellVelZ = cfxCellVelZ[m];
		cfx_cellPressure = cfxCellPressure[m];

		//--------------------------------------------------------
		//cout<<cfdDx<<"  "<<cfdDy<<"  "<<cfdDz<<endl;
		//--Find i,j,k indexs of cell
		int i = (int)floor(cfx_cellX/cfdDx);
		if(i < 0)
		{
        	i = 0;
        }
        else if(i > cfd_iIndex-1)
        {
        	i = cfd_iIndex-1;
        }

		int j = (int)floor(cfx_cellY/cfdDy);
		if(j < 0)
        {
                j = 0;
        }
        else if(j > cfd_jIndex-1)
        {
                j = cfd_jIndex-1;
        }
		int k = (int)floor(cfx_cellZ/cfdDz);
		if(k < 0)
        {
                k = 0;
        }
        else if(k > cfd_kIndex-1)
        {
                k = cfd_kIndex-1;
        }

		//---------------------------

		//-- Now find corresponding array index
		int p = i;
		int q = j*cfd_iIndex;
		int r = k*cfd_iIndex*cfd_jIndex;
		int cfd_cellIndex = p+q+r;

		//if(cellIndex > -1 && cellIndex < iIndex*jIndex*kIndex)
		//{
		//--If the CFD point is closer to cell store the data
		double distance = sqrt(pow((cfdCellCellX[cfd_cellIndex] - (cfx_cellX-cfd_iIndex*0.5*cfdDx)),2)
							+pow((cfdCellCellY[cfd_cellIndex] - (cfx_cellY-cfd_jIndex*0.5*cfdDy)),2)
							+pow((cfdCellCellZ[cfd_cellIndex] - cfx_cellZ),2));
        	//cout<<distance/id->lengthFactor<<"  "<<cfd_neighbourRegion*id->lengthFactor<<endl;
		if (distance < cfd_neighbourRegion*id->lengthFactor)
		{
			cfdCellNeighbourDistance[cfd_cellIndex][cfdCellNoOfPoints[cfd_cellIndex]] = distance;
			cfdCellNeighbourVelX[cfd_cellIndex][cfdCellNoOfPoints[cfd_cellIndex]] = cfx_cellVelX;
			cfdCellNeighbourVelY[cfd_cellIndex][cfdCellNoOfPoints[cfd_cellIndex]] = cfx_cellVelY;
			cfdCellNeighbourVelZ[cfd_cellIndex][cfdCellNoOfPoints[cfd_cellIndex]] = cfx_cellVelZ;
			cfdCellNeighbourPressure[cfd_cellIndex][cfdCellNoOfPoints[cfd_cellIndex]] = cfx_cellPressure;
			cfdCellNoOfPoints[cfd_cellIndex] = cfdCellNoOfPoints[cfd_cellIndex] +1;
			if (cfdCellNoOfPoints[cfd_cellIndex] > neighCfdPoints)
			{
				cout<<"cfdCellNoOfPoints[cfd_cellIndex] > neighCfdPoints"<<endl;
				exit(0);
			}
		}
		//}
		//-------------------------------------------------------------
	}//end of m loop

	cout<<"DEM cells updated"<<endl;

	for (int qq=0; qq<cfd_iIndex*cfd_jIndex*cfd_kIndex; qq++ )
	{
		pre_cfdCellVelX[qq] = cfdCellVelX[qq];
		pre_cfdCellVelY[qq] = cfdCellVelY[qq];
		pre_cfdCellVelZ[qq] = cfdCellVelZ[qq];
		cfdCellVelX[qq] = 0.0;
		cfdCellVelY[qq] = 0.0;
		cfdCellVelZ[qq] = 0.0;
                
		//---- velocity X --------------------------
		double val_0 = cfdCellNeighbourVelX[qq][0];
		double val_1=  cfdCellNeighbourVelX[qq][1];
		double val_2 = cfdCellNeighbourVelX[qq][2];
		double val_3 = cfdCellNeighbourVelX[qq][3];
		double val_4 = cfdCellNeighbourVelX[qq][4];
		
    		double x1 = 0.0;
    		double x2 = 0.0;
    		double x3 = 0.0;
    		double x4 = 0.0;
    		double x5 = 0.0;
             
    		if(cfdCellNeighbourDistance[qq][0] != 0)
    		{
        		x1 = 1.0/cfdCellNeighbourDistance[qq][0];
    		}
    		if(cfdCellNeighbourDistance[qq][1] != 0)
    		{
        		x2 = 1.0/cfdCellNeighbourDistance[qq][1];
    		}
    		if(cfdCellNeighbourDistance[qq][2] != 0)
    		{
        		x3 = 1.0/cfdCellNeighbourDistance[qq][2];
    		}
    		if(cfdCellNeighbourDistance[qq][3] != 0)
    		{
        		x4 = 1.0/cfdCellNeighbourDistance[qq][3];
    		}
    		if(cfdCellNeighbourDistance[qq][4] != 0)
    		{
        		x5 = 1.0/cfdCellNeighbourDistance[qq][4];
    		}
    
		if(x1+x2+x3+x4+x5 != 0)
		{
			//cfdCellVelX[qq] = (1-stepNumber)*cfdCellVelX[qq]/stepNumber + (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5)/stepNumber;
			cfdCellVelX[qq] = (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5);
		}
		//------------------------------------------
		//---- velocity Y --------------------------
		val_0 = cfdCellNeighbourVelY[qq][0];
		val_1=  cfdCellNeighbourVelY[qq][1];
		val_2 = cfdCellNeighbourVelY[qq][2];
		val_3 = cfdCellNeighbourVelY[qq][3];
		val_4 = cfdCellNeighbourVelY[qq][4];

        	if(x1+x2+x3+x4+x5 != 0)
        	{
			//cfdCellVelY[qq] = (1-stepNumber)*cfdCellVelY[qq]/stepNumber + (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5)/stepNumber;
			cfdCellVelY[qq] = (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5);
		}
	
		//------------------------------------------
        //---- velocity Z --------------------------
		val_0 = cfdCellNeighbourVelZ[qq][0];
		val_1=  cfdCellNeighbourVelZ[qq][1];
		val_2 = cfdCellNeighbourVelZ[qq][2];
		val_3 = cfdCellNeighbourVelZ[qq][3];
		val_4 = cfdCellNeighbourVelZ[qq][4];

        	if(x1+x2+x3+x4+x5 != 0)
        	{
			//cfdCellVelZ[qq] = (1-stepNumber)*cfdCellVelZ[qq]/stepNumber + (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5)/stepNumber;
			cfdCellVelZ[qq] = (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5);
		}
        //------------------------------------------
        //---- pressure --------------------------
		val_0 = cfdCellNeighbourPressure[qq][0];
		val_1=  cfdCellNeighbourPressure[qq][1];
		val_2 = cfdCellNeighbourPressure[qq][2];
		val_3 = cfdCellNeighbourPressure[qq][3];
		val_4 = cfdCellNeighbourPressure[qq][4];

        if(x1+x2+x3+x4+x5 != 0)
        {
			//cfdCellPressure[qq] = (1-stepNumber)*cfdCellPressure[qq]/stepNumber + (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5)/stepNumber;
			cfdCellPressure[qq] = (val_0*x1+val_1*x2+val_2*x3+val_3*x4+val_4*x5)/(x1+x2+x3+x4+x5);
        }
        //------------------------------------------
		//---Initialize cellCfdDistance for the next input file -------
        for (int rr=0; rr<neighCfdPoints; rr++)
        {
                cfdCellNeighbourDistance[qq][rr] = 0.0;
                cfdCellNeighbourVelX[qq][rr] = 0.0;
                cfdCellNeighbourVelY[qq][rr] = 0.0;
                cfdCellNeighbourVelZ[qq][rr] = 0.0;
                cfdCellNeighbourPressure[qq][rr] = 0.0;
        }
		
		cfdCellNoOfPoints[qq] = 0; //reset no of neighbour cfx cells
		//----------------------------------------------------------------
	}//end of no of CFD cells
	cfdFile.close();
	
	//----- Find pressure gradient in each cell ----------
	for (int k=0; k<cfd_kIndex; k++)
	{
		for (int j=0; j<cfd_jIndex; j++ )
		{
			for (int i=0; i<cfd_iIndex; i++)
			{
				//------find pressure gradient in X direction-----------------------------
				double pX2 = 0.0;
				double pX1 = 0.0;
				
				int cfdCellIndex = i + (j*cfd_iIndex) + (k*cfd_iIndex*cfd_jIndex);
 		                pre_cfdCellPGradX[cfdCellIndex] = cfdCellPGradX[cfdCellIndex];
                		pre_cfdCellPGradY[cfdCellIndex] = cfdCellPGradY[cfdCellIndex];
                		pre_cfdCellPGradZ[cfdCellIndex] = cfdCellPGradZ[cfdCellIndex];
                		cfdCellPGradX[cfdCellIndex] = 0.0;
                		cfdCellPGradY[cfdCellIndex] = 0.0;
                		cfdCellPGradZ[cfdCellIndex] = 0.0;
 				
				int cfdPostIndex = (i+1) + (j*cfd_iIndex) + (k*cfd_iIndex*cfd_jIndex);
				if ((i+1) > (cfd_iIndex-1))
				{
					pX2 = 0.0;
				}
				else
				{
					pX2 = 0.5*(cfdCellPressure[cfdCellIndex] + cfdCellPressure[cfdPostIndex]);
				}

				int cfdPreIndex = (i-1) + (j*cfd_iIndex) + (k*cfd_iIndex*cfd_jIndex);
				if ((i-1) < 0)
				{
					pX1 = 0.0;
				}
				else
				{
					pX1 = 0.5*(cfdCellPressure[cfdCellIndex] + cfdCellPressure[cfdPreIndex]);
				}
				
				cfdCellPGradX[cfdCellIndex] = (pX2 - pX1)/(id->chamberInnerDia)/(cfd_iIndex*1.0);
				//--------------------------------------------------------------------------

				//--------find pressure gradient in Y direction------------------------------
				double pY2 = 0.0;
				double pY1 = 0.0;
				
				cfdPostIndex = i + (j+1)*cfd_iIndex + (k*cfd_iIndex*cfd_jIndex);
				if ((j+1) > (cfd_jIndex-1))
				{
					pY2 = 0.0;
				}
				else
				{
					pY2 = 0.5*(cfdCellPressure[cfdCellIndex] + cfdCellPressure[cfdPostIndex]);
				}

				cfdPreIndex = i + (j-1)*cfd_iIndex + (k*cfd_iIndex*cfd_jIndex);
				if ((j-1) < 0)
				{
					pY1 = 0.0;
				}
				else
				{
					pY1 = 0.5*(cfdCellPressure[cfdCellIndex] + cfdCellPressure[cfdPreIndex]);
				}
				cfdCellPGradY[cfdCellIndex] = (pY2 - pY1)/(id->chamberInnerDia)/(cfd_jIndex*1.0);
				//-----------------------------------------------------------------------------
				
				//-------find pressure gradient in Z direction--------------------------------
				double pZ2 = 0.0;
				double pZ1 = 0.0;
				
				cfdPostIndex = i + (j*cfd_iIndex) + (k+1)*cfd_iIndex*cfd_jIndex;
				if ((k+1) > (cfd_kIndex-1))
				{
					pZ2 = 0.0;
				}
				else
				{
					pZ2 = 0.5*(cfdCellPressure[cfdCellIndex] + cfdCellPressure[cfdPostIndex]);
				}

				cfdPreIndex = i + (j*cfd_iIndex) + (k-1)*cfd_iIndex*cfd_jIndex;
				if ((k-1) < 0)
				{
					pZ1 = 0.0;
				}
				else
				{
					pZ1 = 0.5*(cfdCellPressure[cfdCellIndex] + cfdCellPressure[cfdPreIndex]);
				}
				cfdCellPGradZ[cfdCellIndex] = (pZ2 - pZ1)/(id->chamberLength)/(cfd_kIndex*1.0);
				//-----------------------------------------------------------------------------
			}
		}
	}

	//----- for testing -----------------------------
	/*ofstream out_File("test_velocity.dat", ios::out);
	if (fcounter == 8)
	{
	out_File<<"ZONE T = 'zone' i="<<cfd_iIndex<<", j= "<<cfd_jIndex<<", k = "<<cfd_kIndex<<endl;
	for (int i=0; i<cfd_iIndex*cfd_jIndex*cfd_kIndex; i++)
	{
		out_File<<setw(15)<<cfdCellCellX[i]<<setw(15)<<cfdCellCellY[i]
			<<setw(15)<<cfdCellCellZ[i]<<setw(15)<<cfdCellPressure[i]
			<<setw(15)<<cfdCellVelX[i]<<setw(15)<<cfdCellVelY[i]<<setw(15)<<cfdCellVelZ[i]<<endl;
	}
	out_File.close();
	exit(0);
	}*/
	//-----------------------------------------------
};
/***************************************************
 * This method finds the contact of the particle with 
 * disk,hole and drum. Call from findNeighbours()
 ***************************************************/
void Mill::findContact()
{	
	for (int i=0;  i<compIndex; i++) //for each component
	{
		findCompContact(i);
	}
};

void Mill::calForceBetweenPart()
{
	if(gap/id->lengthFactor < - gapLimit)
		errorFile<<" par "<<particleNo[parIndex]<<" cellNo  "<<particleCellNo[parIndex]<<" gap "<<gap/id->lengthFactor<<" BETWEEN PARTICLE "<<endl;
	
	ipDia = particleDiameter[parIndex];
	jpDia = particleDiameter[nParIndex];
	rStar = 0.5*ipDia*jpDia/(ipDia+jpDia); 
	temp = particleX[nParIndex];
	ijVec[0] = parX - temp;
	temp = particleY[nParIndex];
	ijVec[1] = parY- temp;
	temp = nParZ;
	ijVec[2] = parZ- temp;

	unitVector(ijVec);
	unitVec[0] = r1;
	unitVec[1] = r2;
	unitVec[2] = r3;

	ipRVec[0] = -0.5*ipDia*unitVec[0];
	ipRVec[1] = -0.5*ipDia*unitVec[1];
	ipRVec[2] = -0.5*ipDia*unitVec[2];

	jpRVec[0] = 0.5*jpDia*unitVec[0];
	jpRVec[1] = 0.5*jpDia*unitVec[1];
	jpRVec[2] = 0.5*jpDia*unitVec[2];
	
	double tempVec1[3];
	tempVec1[0] = particleAngVelX[parIndex];
	tempVec1[1] = particleAngVelY[parIndex];
	tempVec1[2] = particleAngVelZ[parIndex];

	crossProduct(tempVec1,ipRVec);
	rotVel[0] = r1;    
	rotVel[1] = r2;
	rotVel[2] = r3;

	tempVec1[0] = particleVelX[parIndex];
	tempVec1[1] = particleVelY[parIndex];
	tempVec1[2] = particleVelZ[parIndex];

	vecAdd(tempVec1,rotVel);
	ipCntPntVel[0] = r1;    
	ipCntPntVel[1] = r2;
	ipCntPntVel[2] = r3;

	double tempVec2[3];
	tempVec2[0] = particleVelX[nParIndex]; 
	tempVec2[1] = particleVelY[nParIndex]; 
	tempVec2[2] = particleVelZ[nParIndex]; 

	vecSubstract(tempVec1,tempVec2);
	double tempVec3[3];
	tempVec3[0] = r1;
	tempVec3[1] = r2;
	tempVec3[2] = r3;
	nrmVel = dotProduct(tempVec3,unitVec);

	tempVec1[0] = particleAngVelX[nParIndex];
	tempVec1[1] = particleAngVelY[nParIndex];
	tempVec1[2] = particleAngVelZ[nParIndex];

	crossProduct(tempVec1,jpRVec);
	rotVel[0] = r1;    
	rotVel[1] = r2;
	rotVel[2] = r3;
	
	vecAdd(tempVec2,rotVel);
	jpCntPntVel[0] = r1;     
	jpCntPntVel[1] = r2; 
	jpCntPntVel[2] = r3; 

	vecSubstract(ipCntPntVel,jpCntPntVel);
	cntPntVel[0] = r1;    
	cntPntVel[1] = r2;
	cntPntVel[2] = r3;
	instRelativeVel = dotProduct(cntPntVel,unitVec);
	projVec(cntPntVel,unitVec,0); //project relative velocity vector
	tangVel = sqrt(r1*r1+r2*r2+r3*r3);
	sfCoff = id->sfricp;
	normDampC = id->ppNormDampC;
	boundary = 0; //  contact between particles 
	calculateParForce();
};

/*********************************************************
 * Calculate contact forces and momentum for each particle
 * with contact surfaces
 *********************************************************/
void Mill::calculateForce()
{
	totForce[0] = 0.0;
	totForce[1] = 0.0;
	totForce[2] = 0.0;
	momentum[0] = 0.0;
	momentum[1] = 0.0;
	momentum[2] = 0.0;

	nrmDsp = -gap;
	nrmCntForce = id->pElasticMod*sqrt(rStar*nrmDsp)*nrmDsp;
	nrmDampForce = -normDampC*id->pElasticMod*sqrt(rStar*nrmDsp)*nrmVel;//(milliN)scalar//

	tngDampForce = 0;
	contPointDisp[0] = *(cntPntVel)*id->timeStep;
	contPointDisp[1] = *(cntPntVel+1)*id->timeStep;
	contPointDisp[2] = *(cntPntVel+2)*id->timeStep;
	
	projVec(contPointDisp,unitVec,0);// gives the displacement perpendicular to normal vector(shear direction)
	double slidingDisp = sqrt(pow(r1,2)+pow(r2,2)+pow(r3,2)); //sliding distance
	tipContPointDisp[0] = r1;
	tipContPointDisp[1] = r2;
	tipContPointDisp[2] = r3;
	double tempVal = sqrt(pow(r1,2)+pow(r2,2)+pow(r3,2));
	double tempVec1[3];
	tempVec1[0] = 0.0;
	tempVec1[1] = 0.0;
	tempVec1[2] = 0.0;
	if (tempVal != 0)
	{
		tempVec1[0] = r1/tempVal;
		tempVec1[1] = r2/tempVal;
		tempVec1[2] = r3/tempVal;
	}
	else
	{
		tempVec1[0] = 0.0;
		tempVec1[1] = 0.0;
		tempVec1[2] = 0.0;
	}
	double arrayVal[3];
	arrayVal[0] = particleHisDispX[parIndex];
	arrayVal[1] = particleHisDispY[parIndex];
	arrayVal[2] = particleHisDispZ[parIndex];

	projVec(arrayVal,unitVec,1);
	disp[0] = r1;
	disp[1] = r2; 
	disp[2] = r3; 
	disp[0] = *disp + *tipContPointDisp;
	disp[1] = *(disp+1) + *(tipContPointDisp+1);
	disp[2] = *(disp+2) + *(tipContPointDisp+2);
	dsmax = dsmaxCff*nrmDsp;

	dd = sqrt(pow(*disp,2)+pow(*(disp+1),2)+pow(*(disp+2),2));
	fdt = sfCoff*nrmCntForce;

	if (dd < 1e-6)
	{
		dd = 0.0;
	}
	if (dti < 1e-6)
	{
		dti = 0.0;
	}

	if (dd >= dsmax)
	{
		disp[0] = *(disp)*dsmax/dd;
		disp[1] = *(disp+1)*dsmax/dd;
		disp[2] = *(disp+2)*dsmax/dd;
		dti = sqrt(pow(*tipContPointDisp,2)+pow(*(tipContPointDisp+1),2)+pow(*(tipContPointDisp+2),2));
		if (dti != 0)
		{
			fdtx = - fdt*(*tipContPointDisp)/dti;
			fdty = - fdt*(*(tipContPointDisp+1))/dti;
			fdtz = - fdt*(*(tipContPointDisp+2))/dti;
		}
		else//if dti = 0, then tipContPointDisp also zero. Therefor ftdx,fty,ftz should be zero.
		{
			fdtx = 0.0;
			fdty = 0.0;
			fdtz = 0.0;
		}
	}
	else																																																																																																																																																																					
	{
		if (dd != 0)
		{
			fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
			fdtx = -fdt*disp[0]/dd;
			fdty = -fdt*disp[1]/dd;
			fdtz = -fdt*disp[2]/dd;
		}
		else
		{
			fdtx = 0.0;
			fdty = 0.0;
			fdtz = 0.0;
		}
	}
	
	particleHisDispX[parIndex] = disp[0];
	particleHisDispY[parIndex] = disp[1];
	particleHisDispZ[parIndex] = disp[2];

	//sum of forces
	nrmForce = (nrmCntForce + nrmDampForce);
	temp = nrmForce*(*unitVec);
	totForce[0] = fdtx + temp;
	temp = nrmForce*(*(unitVec+1));
	totForce[1] = fdty + temp;
	temp = nrmForce*(*(unitVec+2));
	totForce[2] = fdtz + temp;

	projVec(totForce,unitVec,0); // gives the projected force vector on a plane perpendicular to normal vector
	instShearForce = sqrt(pow(r1,2) + pow(r2,2) + pow(r3,2)); 

	instCompForce = totForce[0]*unitVec[0]+totForce[1]*unitVec[1]+totForce[2]*unitVec[2];
	//particleInstantCompForce[parIndex] = instCompForce;
	//particleInstantShearForce[parIndex] = instShearForce;
	//particleAverageCompForce[parIndex] = particleAverageCompForce[parIndex] + instCompForce;
	//particleAverageShearForce[parIndex] = particleAverageShearForce[parIndex] + instShearForce;
	//particleAverageCollE[parIndex] =  particleAverageCollE[parIndex] + collEnergy;
    particleNoOfContacts[parIndex] = particleNoOfContacts[parIndex] + 1;
	instImpactEnergy = (-nrmDampForce*nrmVel*id->timeStep+sqrt(fdtx*fdtx+fdty*fdty+fdtz*fdtz)*slidingDisp)/id->energyFactor;
    instShearEnergy = instShearForce*slidingDisp/id->energyFactor;
	collEnergy =  0.5*(particleMass[parIndex]*instRelativeVel*instRelativeVel)/id->energyFactor;
	tangCollEnergy = 0.5*(particleMass[parIndex]*tangVel*tangVel)/id->energyFactor;
	collPointX = particleX[parIndex]-unitVec[0]*particleDiameter[parIndex]*0.5; // collision point
	collPointY = particleY[parIndex]-unitVec[1]*particleDiameter[parIndex]*0.5;
	collPointZ = particleZ[parIndex]-unitVec[2]*particleDiameter[parIndex]*0.5;
   
	double wear1 = 0.0;
/*	
	//------- Wear calculation ------------------------------------------------------
    if (collPointZ < zMaximum*id->lengthFactor && collPointZ > zMinimum*id->lengthFactor)
    {
		double wearAng = fabs(M_PI*0.5 - acos((cntPntVel[0]*unitVec[0]+cntPntVel[1]*unitVec[1]+cntPntVel[2]*unitVec[2])
								/(sqrt(pow(cntPntVel[0],2)+pow(cntPntVel[1],2)+pow(cntPntVel[2],2))*sqrt(pow(unitVec[0],2)+pow(unitVec[1],2)+pow(unitVec[2],2)))));

		if(wearAng*180.0/M_PI > 18.5)
        {
            wear1 = collEnergy*pow(cos(wearAng),2)/2.0;
        }
        else
        {
            wear1 =  collEnergy*(sin(2*wearAng) - 3*pow(sin(wearAng),2))/2.0;
        }

    }
	//----------- End of wear calculation -------------------------------------------

*/
    if (checkData(1, parIndex, boundary, collPointX,collPointY,collPointZ,
						 instCompForce, instImpactEnergy,instShearEnergy, gap, collEnergy,tangCollEnergy,wear1) == false)

	{
		//increase no of collisions by one
		//cellArrayList[particleCellNo[parIndex] - 1]->increaseCollNo();
		totalWear1 = totalWear1 + wear1;
	}
	totalWear2 = totalWear2 + instShearEnergy/id->energyFactor;
	/*  
	if (collPointZ < zMaximum*id->lengthFactor && collPointZ > zMinimum*id->lengthFactor)
	{
		for (int i=0; i<diskFileArraySize1; i++)
        {
			if (sqrt(pow((collPointX - pointArrayX[i]),2)+pow((collPointY - pointArrayY[i]),2)+pow((collPointZ - pointArrayZ[i]),2)) < wearGapLimit*id->lengthFactor)
			{
				pointData1[i] = pointData1[i] + totalWear[parIndex];
				pointData2[i] = pointData2[i] + instShearEnergy/id->energyFactor;
			}
		}
	}*/
    totalWear[parIndex] = 0.0; //reset wear data

	collPointX = particleX[parIndex];// X coordinate of the particle
	collPointY = particleY[parIndex]; // Y coordinate of the particle
	collPointZ = particleZ[parIndex]; // Z coordinate of the particle
	        
    rotorForce[0] = -nrmCntForce*(*unitVec);
	rotorForce[1] = -nrmCntForce*(*(unitVec+1));
	rotorForce[2] = -nrmCntForce*(*(unitVec+2));
		
	crossProduct(ipRVec,totForce);
	momentum[0] = r1;
	momentum[1] = r2;
	momentum[2] = r3;
	double velArray[3];
	velArray[0] = particleAngVelX[parIndex];
	velArray[1] = particleAngVelY[parIndex];
	velArray[2] = particleAngVelZ[parIndex];
	momentumR[0] = (*velArray)*(-id->rollfric*parDia*nrmCntForce);
	momentumR[1] = (*(velArray+1))*(-id->rollfric*parDia*nrmCntForce);
	momentumR[2] = (*(velArray+2))*(-id->rollfric*parDia*nrmCntForce);
	vecAdd(momentumR,momentum);
	momentum[0] = r1;     
	momentum[1] = r2; 
	momentum[2] = r3; 
	
	vecAdd(totForce,pForce);//total force on the particle
	pForce[0] = r1;
	pForce[1] = r2;
	pForce[2] = r3;
	vecAdd(pMom,momentum);
	pMom[0] = r1;   //total momentum on the particle
	pMom[1] = r2; 
	pMom[2] = r3;
	
	particleForceX[parIndex] = pForce[0];
	particleForceY[parIndex] = pForce[1];
	particleForceZ[parIndex] = pForce[2];
	particleMomentumX[parIndex] = pMom[0];
	particleMomentumY[parIndex] = pMom[1];
	particleMomentumZ[parIndex] = pMom[2];
};//end of calculateForce()


/************************************************************
 * Calculates particle-particle contact forces and momentum 
 *************************************************************/
void Mill::calculateParForce()
{
	totForce[0] = 0.0;
	totForce[1] = 0.0;
	totForce[2] = 0.0;
	momentum[0] = 0.0;
	momentum[1] = 0.0;
	momentum[2] = 0.0;
	
	nrmDsp = -gap;
	nrmCntForce = id->pElasticMod*sqrt(rStar*nrmDsp)*nrmDsp;
	nrmDampForce = -normDampC*id->pElasticMod*sqrt(rStar*nrmDsp)*nrmVel;//(milliN)scalar//

	tngDampForce = 0;
	
	contPointDisp[0] = *(cntPntVel)*id->timeStep;
	contPointDisp[1] = *(cntPntVel+1)*id->timeStep;
	contPointDisp[2] = *(cntPntVel+2)*id->timeStep;
	
	projVec(contPointDisp,unitVec,0);// gives the displacement perpendicular to normal vector(shear direction)
	double slidingDisp = sqrt(pow(r1,2)+pow(r2,2)+pow(r3,2)); //sliding distance
	tipContPointDisp[0] = r1;
	tipContPointDisp[1] = r2;
	tipContPointDisp[2] = r3;
	double tempVal1 = sqrt(pow(r1,2)+pow(r2,2)+pow(r3,2));
	double tempVec1[3];
	tempVec1[0] = 0.0;
	tempVec1[1] = 0.0;
	tempVec1[2] = 0.0;
	if (tempVal != 0)
	{
		tempVec1[0] = r1/tempVal1;
		tempVec1[1] = r2/tempVal1;
		tempVec1[2] = r3/tempVal1;
	}
	else
	{
		tempVec1[0] = 0.0;
		tempVec1[1] = 0.0;
		tempVec1[2] = 0.0;
	}
	double arrayVal[3];
	arrayVal[0] = particleHisDispX[parIndex];
	arrayVal[1] = particleHisDispY[parIndex];
	arrayVal[2] = particleHisDispZ[parIndex];

	projVec(arrayVal,unitVec,1);
	disp[0] = r1;
	disp[1] = r2; 
	disp[2] = r3; 
	disp[0] = *disp + *tipContPointDisp;
	disp[1] = *(disp+1) + *(tipContPointDisp+1);
	disp[2] = *(disp+2) + *(tipContPointDisp+2);

	dsmax = dsmaxCff*nrmDsp;

	dd = sqrt(pow(*disp,2)+pow(*(disp+1),2)+pow(*(disp+2),2));
	fdt = sfCoff*nrmCntForce;

	if (dd < 1e-6)
	{
		dd = 0.0;
	}
	if (dti < 1e-6)
	{
		dti = 0.0;
	}

	if (dd >= dsmax)
	{
		disp[0] = *(disp)*dsmax/dd;
		disp[1] = *(disp+1)*dsmax/dd;
		disp[2] = *(disp+2)*dsmax/dd;
		dti = sqrt(pow(*tipContPointDisp,2)+pow(*(tipContPointDisp+1),2)+pow(*(tipContPointDisp+2),2));
		if (dti != 0)
		{
			fdtx = - fdt*(*tipContPointDisp)/dti;
			fdty = - fdt*(*(tipContPointDisp+1))/dti;
			fdtz = - fdt*(*(tipContPointDisp+2))/dti;
		}
		else//if dti = 0, then tipContPointDisp also zero. Therefor ftdx,fty,ftz should be zero.
		{
			fdtx = 0.0;
			fdty = 0.0;
			fdtz = 0.0;
		}
	}
	else																																																																																																																																																																					
	{
		if (dd != 0)
		{
			fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
			fdtx = -fdt*disp[0]/dd;
			fdty = -fdt*disp[1]/dd;
			fdtz = -fdt*disp[2]/dd;
			//File<<" dd != 0 dti "<<dti<<" fdtx "<<fdtx<<" fdty "<<fdty<<" fdtz "<<fdtz;
		}
		else
		{
			fdtx = 0.0;
			fdty = 0.0;
			fdtz = 0.0;
		}
	}

	particleHisDispX[parIndex] = disp[0];
	particleHisDispY[parIndex] = disp[1];
	particleHisDispZ[parIndex] = disp[2];

	//sum of forces
	nrmForce = (nrmCntForce + nrmDampForce);
	temp = nrmForce*(*unitVec);
	totForce[0] = fdtx + temp;
	temp = nrmForce*(*(unitVec+1));
	totForce[1] = fdty + temp;
	temp = nrmForce*(*(unitVec+2));
	totForce[2] = fdtz + temp;
	// Sum of force components in tangential direction //
	projVec(totForce,unitVec,0); // gives the projected force vector on a plane perpendicular to normal vector
	instShearForce = sqrt(pow(r1,2) + pow(r2,2) + pow(r3,2)); // (/1e-6 N)
	// find a unit vector along the projected vector
	// Now the resultant of the project vector gives the compressive force
	instCompForce = totForce[0]*unitVec[0]+totForce[1]*unitVec[1]+totForce[2]*unitVec[2];
	instShearEnergy = 0.5*instShearForce*slidingDisp;
	//par->setInstantForces(instCompForce,instShearForce);
	//particleInstantCompForce[parIndex] = instCompForce;
	//particleInstantShearForce[parIndex] = instShearForce;
	//particleInstantCompForce[nParIndex] = instCompForce;
	//particleInstantShearForce[nParIndex] = instShearForce;
	//particleAverageCompForce[parIndex] = particleAverageCompForce[parIndex] + instCompForce;
	//particleAverageShearForce[parIndex] = particleAverageShearForce[parIndex] + instShearForce;
	// particleAverageCollE[parIndex] =  particleAverageCollE[parIndex] + collEnergy;
    particleNoOfContacts[parIndex] = particleNoOfContacts[parIndex] + 1;
	//particleAverageCompForce[nParIndex] = particleAverageCompForce[nParIndex] + instCompForce;
	//particleAverageShearForce[nParIndex] = particleAverageShearForce[nParIndex] + instShearForce;
	//particleAverageCollE[nParIndex] =  particleAverageCollE[nParIndex] + collEnergy;
    particleNoOfContacts[nParIndex] = particleNoOfContacts[nParIndex] + 1;

	rotorForce[0] = -nrmCntForce*(*unitVec);
	rotorForce[1] = -nrmCntForce*(*(unitVec+1));
	rotorForce[2] = -nrmCntForce*(*(unitVec+2));
		
	crossProduct(ipRVec,totForce);
	momentum[0] = r1;
	momentum[1] = r2;
	momentum[2] = r3;
	
	double velArray[3];
	velArray[0] = particleAngVelX[parIndex];
	velArray[1] = particleAngVelY[parIndex];
	velArray[2] = particleAngVelZ[parIndex];

	momentumR[0] = (*velArray)*(-id->rollfric*parDia*nrmCntForce);
	momentumR[1] = (*(velArray+1))*(-id->rollfric*parDia*nrmCntForce);
	momentumR[2] = (*(velArray+2))*(-id->rollfric*parDia*nrmCntForce);
	vecAdd(momentumR,momentum);
	momentum[0] = r1;     
	momentum[1] = r2; 
	momentum[2] = r3; 
	
	vecAdd(totForce,pForce);//total force on the particle
	pForce[0] = r1;
	pForce[1] = r2;
	pForce[2] = r3;
	vecAdd(pMom,momentum);
	pMom[0] = r1;   //total momentum on the particle
	pMom[1] = r2; 
	pMom[2] = r3;
	
	particleForceX[parIndex] = pForce[0];
	particleForceY[parIndex] = pForce[1];
	particleForceZ[parIndex] = pForce[2];
	particleMomentumX[parIndex] = pMom[0];
	particleMomentumY[parIndex] = pMom[1];
	particleMomentumZ[parIndex] = pMom[2];
	instImpactEnergy = (-nrmDampForce*nrmVel*id->timeStep+sqrt(fdtx*fdtx+fdty*fdty+fdtz*fdtz)*slidingDisp)/id->energyFactor;
    instShearEnergy = instShearForce*slidingDisp/id->energyFactor;
	collEnergy =  0.5*(0.5*(particleMass[parIndex]+particleMass[nParIndex])*instRelativeVel*instRelativeVel)/id->energyFactor;
	tangCollEnergy = 0.5*(0.5*(particleMass[parIndex]+particleMass[nParIndex])*tangVel*tangVel)/id->energyFactor;

	collPointX = particleX[parIndex]+ipRVec[0]; //X coordinate of the collision point
	collPointY = particleY[parIndex]+ipRVec[1]; //Y coordinate of the collision point
	collPointZ = particleZ[parIndex]+ipRVec[2]; //Z coordinate of the collision point
    checkData(1, parIndex, particleArrayIndex[nParIndex], collPointX,collPointY,collPointZ,
		instCompForce, instImpactEnergy,instShearEnergy, gap, collEnergy,tangCollEnergy, 0.0);
	
	//particleData[parIndex]->checkData(particleArrayIndex[nParIndex], collPointX,collPointY,collPointZ,
	//	instCompForce, instShearForce,instShearEnergy, gap, collEnergy,instRelativeVel);
	collPointX = particleX[parIndex];// X coordinate of the particle
	collPointY = particleY[parIndex]; // Y coordinate of the particle
	collPointZ = particleZ[parIndex]; // Z coordinate of the particle

	// Force and Momentum on the other particle //
	totForce[0] = *(totForce)*(-1.0);
	totForce[1] = *(totForce+1)*(-1.0);
	totForce[2] = *(totForce+2)*(-1.0);

	particleForceX[nParIndex] = particleForceX[nParIndex] + totForce[0];
	particleForceY[nParIndex] = particleForceY[nParIndex] + totForce[1];
	particleForceZ[nParIndex] = particleForceZ[nParIndex] + totForce[2];

    crossProduct(jpRVec,totForce);
    momentum[0] = r1;
    momentum[1] = r2;
    momentum[2] = r3;

    double nVelArray[3];
    nVelArray[0] = particleAngVelX[nParIndex];
    nVelArray[1] = particleAngVelY[nParIndex];
    nVelArray[2] = particleAngVelZ[nParIndex];

    momentumR[0] = (*nVelArray)*(-id->rollfric*particleDiameter[nParIndex]*nrmCntForce);
    momentumR[1] = (*(nVelArray+1))*(-id->rollfric*particleDiameter[nParIndex]*nrmCntForce);
    momentumR[2] = (*(nVelArray+2))*(-id->rollfric*particleDiameter[nParIndex]*nrmCntForce);
    vecAdd(momentumR,momentum);
    momentum[0] = r1;
    momentum[1] = r2;
    momentum[2] = r3;
	particleMomentumX[nParIndex] = particleMomentumX[nParIndex] + momentum[0];
	particleMomentumY[nParIndex] = particleMomentumY[nParIndex] + momentum[1];
	particleMomentumZ[nParIndex] = particleMomentumZ[nParIndex] + momentum[2];
	checkData(0, nParIndex, particleArrayIndex[parIndex], collPointX,collPointY,collPointZ,
		instCompForce, instImpactEnergy,instShearEnergy, gap, collEnergy,tangCollEnergy,0.0);

	unitVec[0] = *(unitVec)*(-1.0);
	unitVec[1] = *(unitVec+1)*(-1.0);
	unitVec[2] = *(unitVec+2)*(-1.0);

	contPointDisp[0] = *(contPointDisp)*(-1.0);
	contPointDisp[1] = *(contPointDisp+1)*(-1.0);
	contPointDisp[2] = *(contPointDisp+2)*(-1.0);

	projVec(contPointDisp,unitVec,0);

	tipContPointDisp[0] = r1;
	tipContPointDisp[1] = r2;
	tipContPointDisp[2] = r3;
	
	double arrayVal1[3];
	arrayVal1[0] = particleHisDispX[nParIndex];
	arrayVal1[1] = particleHisDispY[nParIndex];
	arrayVal1[2] = particleHisDispZ[nParIndex];
	projVec(arrayVal1,unitVec,1);
	disp[0] = r1;
	disp[1] = r2; 
	disp[2] = r3; 
	disp[0] = *disp + *tipContPointDisp;
	disp[1] = *(disp+1) + *(tipContPointDisp+1);
	disp[2] = *(disp+2) + *(tipContPointDisp+2);

	dd = sqrt(pow(*disp,2)+pow(*(disp+1),2)+pow(*(disp+2),2));
	
	if (dd < 1e-6)
	{
		dd = 0.0;
	}

	if (dd >= dsmax)
	{
		disp[0] = *(disp)*dsmax/dd;
		disp[1] = *(disp+1)*dsmax/dd;
		disp[2] = *(disp+2)*dsmax/dd;
	}
	particleHisDispX[nParIndex] = disp[0];
	particleHisDispY[nParIndex] = disp[1];
	particleHisDispZ[nParIndex] = disp[2];

};//end of calculateParForce()

int Mill::vecAdd(double *v1,double *v2)
{	
	r1 = *v1+ *v2;
	r2 = *(v1+1)+ *(v2+1);
	r3 = *(v1+2)+ *(v2+2);
	return 1;
};


int Mill::vecSubstract(double *v1,double *v2)
{
	r1 = *v1 - *v2;
	r2 = *(v1+1) - *(v2+1);
	r3 = *(v1+2) - *(v2+2);
	return 1;
};


int Mill::crossProduct(double *v1,double *v2)
{
	val1 = *(v1+1)*(*(v2+2));
	val2 = *(v2+1)*(*(v1+2));
	r1 = val1 - val2;
	val1 = *v2*(*(v1+2));
	val2 = *v1*(*(v2+2));
	r2 = val1 - val2;
	val1 = *v1*(*(v2+1));
	val2 = *v2*(*(v1+1));
	r3 = val1 - val2;
	return 1;
};

int Mill::unitVector(double *v)
{
	if(vecMag(v) > 0)
	{
		r1 = *v/sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
		r2 = *(v+1)/sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
		r3 = *(v+2)/sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
	}
	else
	{
		r1 = 0.0;
		r2 = 0.0;
		r3 = 0.0;
	}
	return 1;

};

double Mill::dotProduct(double *v1,double *v2)
{
	temp = *v1*(*v2) + *(v1+1)*(*(v2+1)) + *(v1+2)*(*(v2+2));
	return temp;
};

double Mill::vecMag(double *v)
{
	temp = sqrt(pow(*v,2)+pow(*(v+1),2)+pow(*(v+2),2));
	return temp;
};

int Mill::sclVecMult(double scl, double *v)
{
	r1 = scl*(*v);
	r2 = scl*(*(v+1));
	r3 = scl*(*(v+2));
	return 1;
};

int Mill::projVec(double *v,double *n,int type)
{	if (type == 0)
	{
		crossProduct(n,v);
		tempVec[0] = r1;
		tempVec[1] = r2;
		tempVec[2] = r3;
		
		vec[0] = -(*n);
		vec[1] = -(*(n+1));
		vec[2] = -(*(n+2));
		crossProduct(vec,tempVec);
		tempVec[0] = r1;
		tempVec[1] = r2;
		tempVec[2] = r3;
		return 1;
	}
	else
	{
		crossProduct(n,v);
		tempVec[0] = r1;
		tempVec[1] = r2;
		tempVec[2] = r3;
		vec[0] = -(*n);     
		vec[1] = -(*(n+1)); 
		vec[2] = -(*(n+2)); 
		crossProduct(vec,tempVec);
		tempVec[0] = r1;
		tempVec[1] = r2;
		tempVec[2] = r3;
		temp = pow(*v,2)+pow(*(v+1),2)+pow(*(v+2),2);
		/*
		if (temp < 0.001)
		{
			return 1;
		}*/
		tempVal = sqrt((pow(*tempVec,2)+pow(*(tempVec+1),2)+pow(*(tempVec+2),2))
										/(pow(*v,2)+pow(*(v+1),2)+pow(*(v+2),2)));
		tempVec[0] = (*tempVec)*tempVal;
		tempVec[1] = (*(tempVec+1))*tempVal;
		tempVec[2] = (*(tempVec+2))*tempVal;
		return 1;
	}
};

