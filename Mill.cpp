/**************************************************************
 * Mill.cpp
 * Driving class of the program.
 * Copyright: Copyright (c) 2007
 * SIMPAS - University of New South Wales, Sydney
 * @author Chandana Jayasundara 		
 * @version 5.0
 * 
 **************************************************************/
#include "Mill.h"
#define PI 3.14159

int main()
{
	Mill *mill = new Mill();
	mill->run();
	delete mill;
	system("PAUSE");
	return 0;
};

Mill::Mill()
{
	stepsBetweenFrames = 0;
	firstStep = true;
	ch = '"'; 
	stepNumber = 1;
	noOfCycles = 0;
	frameStepCounter = 0;
	gapLimit = 0.0;
	maxTime = totalTime = startingTime = 0.0;
	collisionGapCoeff = 0.0;
	firstInsertion = true;
	firstTime = true;
	insertable = false;
	counter = tempColor = count = fileCounter = arrCounter = timeCounter = 0;
	zoneCounter = 1;
	cellParIndex = 0;
	cellNo = 0;
	tempX = tempY = tempZ = tempR = dxDot = dyDot = dzDot = xAngVel = yAngVel = zAngVel = 0.0;
	dist = 0.0;
	rotorSpeed = tempTheta = 0.0;
	noOfParticles = index = 0;
	rotorAngPosition = 0.0;
	tempAngPosition = 0.0;
	diskAngPosition = 0.0;
	
	id =  NULL;
	ddR = ddTheta = ddZ = 0.0;
	noOfParticles = 0;
	xCoordinate = yCoordinate = zCoordinate = 0.0;
	dx = dy = dz = cfdDx = cfdDy = cfdDz = pX = pY = pZ = cellDx = cellDy = cellDz = 0.0;
	diskY = 0.0;

	dsmaxCff = tempDist = 0.0;
	contact = -1;
	timeCount = timeCount1 = 0;
	torque = totalTorque = 0.0;
    totTime = 0.0;
	//energyLoss = strainEnergy = 0.0;
	nrmForce = nrmDsp = tngDampForce = instRelativeVel = 0.0;
	dsmax = dd = fdt = dti = fdtx = fdty = fdtz = 0.0;   
	instShearForce = instImpactEnergy = instCompForce = instShearEnergy = collEnergy = tangCollEnergy = collPointX = collPointY = collPointZ  = 0.0;

	gForce = diskThick = diskDia = shaftDia = 0.0;
	pXY = 0.0;
	sfCoff = normDampC = 0.0;
	r1 = r2 = r3 = ipDia = jpDia = rStar = nrmVel = tangVel = 0.0;
	ipX = ipY = ipZ = xx = yy = zz = xy = 0.0;
	jpX = jpY = jpZ = nrmCntForce = nrmDampForce = 0.0;

	for (int i=0; i<3 ;i++ )
	{
		tempVec[i] = 0.0;
		vec[i] = 0.0;
		tipContPointDisp[i] = 0.0;
		disp[i] = 0.0;          
		contPointDisp[i] = 0.0; 
		unitVec[i] = 0.0;       
		ipRVec[i] = 0.0;        
		jpRVec[i] = 0.0;        
		ijVec[i] = 0.0;         
		rotVel[i] = 0.0;        
		ipCntPntVel[i] = 0.0;   
		jpCntPntVel[i] = 0.0;   
		cntPntVel[i] = 0.0;     
		totForce[i] = 0.0;      
		rotorForce[i] = 0.0;    
		rotorTotForce[i] = 0.0; 
		pMom[i] = 0.0;          
		torqueOnMill[i] = 0.0;  
		momentum[i] = 0.0;      
		momentumR[i] = 0.0;     
		pForce[i] = 0.0;        
		drumRVec[i] = 0.0;      
		shaftRVec[i] = 0.0;     
		diskRVec[i] = 0.0;      
		holeRVec[i] = 0.0;      
		shaftCntPntVel[i] = 0.0;
		diskCntPntVel[i] = 0.0; 
		shaftAngVel[i] = 0.0; 
	}

	for (int i=0; i<compArraySize ;i++ )
	{
		componentX[i] = 0.0;
		componentY[i] = 0.0;
		componentZ[i] = 0.0;
		componentR[i] = 0.0;
		componentDia[i] = 0.0;
		componentTheta[i] = 0.0;
		componentLength[i] = 0.0;
		componentHeight[i] = 0.0;
		componentWidth[i] = 0.0;
		noOfSubComponents[i] = 0;   
		componentType[i] = 0;
		componentAngVelocity[i] = 0.0;
		componentSlidingFriction[i] = 0.0;
		for (int j=0; j<hollowCompArraySize ;j++ )
		{
			hollowCompR[i][j] = 0.0;     
			hollowCompTheta[i][j] = 0.0; 
			hollowCompZ[i][j] = 0.0;     
			hollowCompLength[i][j] = 0.0;
			hollowCompDia[i][j] = 0.0;  
			hollowCompAngVelocity[i][j] = 0.0;
		}
	}
	isInside = false;  
	calculateTorque = false;
	isContacted = 0;
	nearestSpokeIndex = 0;
	compIndex = 0;
	parXDash = parYDash = parZDash = 0.0;
	nearestHoleX = nearestHoleY = nearestDiskZ = nearestHubZ = 0.0;
	shaftIndex = drumIndex = firstEndWallIndex = secondEndWallIndex = 0;
	nParIndex = parIndex =  num = boundary = 0;
	parX = parY = parZ = parDia = gap = collisionGap = nParZ = 0.0;
	val1 = val2 = val3 = temp = tempVal = 0.0;
	dragFX = dragFY = dragFZ = bouyancyF = 0.0;
	compType = 0; 
	hubDia = hubThick = 0.0;
	torqueRadius = 0.0;

	for (int i=0; i<particleArraySize; i++ )
	{
		particleNo[i] = 0;
		particleCellNo[i] = 0;
		particleArrayIndex[i] = 0;
		particleDiameter[i] = 0.0;
		particleDensity[i] = 0.0;
		particleColor[i] = 0;
		particleNoOfContacts[i]= 0;
		particleInsertable[i] = true;
		particleContactCalculated[i] = false;
		particleX[i] = 0.0;
		particleY[i] = 0.0;
		particleZ[i] = 0.0;

		particleVelX[i] = 0.0;
		particleVelY[i] = 0.0;
		particleVelZ[i] = 0.0;
		
		particleHisDispX[i] = 0.0;
		particleHisDispY[i] = 0.0;
		particleHisDispZ[i] = 0.0;
		particleDisplacement[i] = 0.0;
		particleAngVelX[i] = 0.0;
		particleAngVelY[i] = 0.0;
		particleAngVelZ[i] = 0.0;
		//particleAverageCompForce[i] = 0.0;
		//particleAverageShearForce[i] = 0.0;
		//particleAverageCollE[i] = 0.0;
        //particleInstantCompForce[i] = 0.0;
		//particleInstantShearForce[i] = 0.0;
		particleMass[i] = 0.0;
		particleInertia[i] = 0.0;
		particleForceX[i] = 0.0;
		particleForceY[i] = 0.0;
		particleForceZ[i] = 0.0;
	    	
		particleDragForceX[i] = 0.0;
        	particleDragForceY[i] = 0.0;
        	particleDragForceZ[i] = 0.0;
                particlePForceX[i] = 0.0;
                particlePForceY[i] = 0.0;
                particlePForceZ[i] = 0.0;
		particleReynolds[i] = 0.0;

		particleMomentumX[i] = 0.0;
		particleMomentumY[i] = 0.0;
		particleMomentumZ[i] = 0.0;
	
		noOfContacts[i] = 0;                        
		totalCollisions[i] = 0;   
		neighListSize[i] = 0;
	
		instantContactListSize[i] = 0;
	
		for (int j=0; j<collisionIndex; j++ )
		{
			instantContactList[i][j] = 0; 
		}

		for (int j=0; j<neighListArraySz; j++ )
		{
			neighbourList[i][j] = 0;     
		}
		for (int j=0; j<arraySz; j++ )
		{
			contactList[i][j] = 0;     
			maxKineticEnergy[i][j] = 0.0;
			collisionEnergy[i][j] = 0.0;
			tangCollisoinEnergy[i][j] = 0.0;
			impactEnergyDissipation[i][j] = 0.0;
			shearEnergyDissipation[i][j] = 0.0;
			previousOverlap[i][j] = 0.0;
		}

		DRUM_OUTER_FACE[i] = false;
		DRUM_SIDE_FACE[i] = false;
		SHAFT[i] = false;
		DISK_SIDE_FACE[i] = false;
		HUB_OUTER_FACE[i] = false;
		HUB_SIDE_FACE[i] = false;	

		totalWear[i] = 0.0;
	}
	
	for (int i=0; i<noOfHoleDivisions ;i++ )
	{
		holeWearData1[i] = 0.0;
		holeWearData2[i] = 0.0;
	}
	totalWear1 = 0.0;
	totalWear2 = 0.0;
	
	cellVolume = 0.0;
	for (int i=0; i<iIndex*jIndex*kIndex; i++)
	{
		cellCellNo[i] = 0;                                                    
		cellNoOfContacts[i] = 0;
		cellCellX[i] = 0.0;                                                     
		cellCellY[i] = 0.0;                                                     
		cellCellZ[i] = 0.0;                                                     
		cellImpactEnergy[i] = 0.0;                                                      
		cellShearEnergy[i] = 0.0;                                                      
		cellCollisionEnergy[i] = 0.0; 
		cellTangCollisionEnergy[i] = 0.0; 
		cellPorosity[i] = 0.0;   
		cellCellListSize[i] = 0;
		cellParticleListSize[i] = 0;	
		
		for (int j=0; j<neighCellListSize; j++)
		{
			cellCellList[i][j] = 0;
		}

		for (int j=0; j<neighParListSize; j++)
		{
			cellParticleList[i][j] = 0;
		}
	}

	for (int i=0; i<cfd_iIndex*cfd_jIndex*cfd_kIndex ;i++ )
	{
		cfdCellNoOfPoints[i] = 0;
		cfdCellCellX[i] = 0.0;
		cfdCellCellY[i] = 0.0;
		cfdCellCellZ[i] = 0.0;
						
		cfdCellVelX[i] = 0.0;
		cfdCellVelY[i] = 0.0;
		cfdCellVelZ[i] = 0.0;
		cfdCellPGradX[i] = 0.0;
		cfdCellPGradY[i] = 0.0;
		cfdCellPGradZ[i] = 0.0;
		cfdCellPressure[i] = 0.0;
		pre_cfdCellVelX[i] = 0.0;
                pre_cfdCellVelY[i] = 0.0;
                pre_cfdCellVelZ[i] = 0.0;
                pre_cfdCellPGradX[i] = 0.0;
                pre_cfdCellPGradY[i] = 0.0;
                pre_cfdCellPGradZ[i] = 0.0;
		
		for (int j=0; j<neighCfdPoints;j++ )
		{
			cfdCellNeighbourDistance[i][j] = 0.0;
			cfdCellNeighbourPressure[i][j] = 0.0;
			cfdCellNeighbourVelX[i][j] = 0.0;
			cfdCellNeighbourVelY[i][j] = 0.0;
			cfdCellNeighbourVelZ[i][j] = 0.0;
		}
	}
	for (int i=0; i<cfxCellPoints; i++ )
	{
		cfxCellX[i] = 0.0;
		cfxCellY[i] = 0.0;
		cfxCellZ[i] = 0.0;
		tempCfxCellX[i] = 0.0;
		tempCfxCellY[i] = 0.0;
		cfxCellVelX[i] = 0.0;
		cfxCellVelY[i] = 0.0;
		cfxCellVelZ[i] = 0.0;
		cfxCellPressure[i] = 0.0;

	}

	cfdFrameCounter = 0;
	cfdFrameAngle = 0;

	/*for (int i=0; i<cfxCellPoints;i++ )
	{
		cfxCellX[i] = 0.0;
		cfxCellY[i] = 0.0;
		cfxCellZ[i] = 0.0;
		cfxCellVelX[i] = 0.0;
		cfxCellVelY[i] = 0.0;
		cfxCellVelZ[i] = 0.0;
		cfxCellPressure[i] = 0.0;
	}*/

	for (int i=0; i<27 ;i++ )
	{
		neighCellList[i] = 0;
	}

	for (int i=0; i<cellParListSize ;i++ )
	{	
		parList[i] = 0;
	}

	for (int i=0; i<diskFileArraySize1 ;i++ )
	{
		pointArrayX[i] = 0.0; 
		pointArrayY[i] = 0.0; 
		pointArrayZ[i] = 0.0; 
		pointData1[i] = 0.0;
		pointData2[i] = 0.0;
	}

	for (int i=0; i<diskFileArraySize2 ;i++ )
	{
		elementArrayX[i] = 0.0; 
		elementArrayY[i] = 0.0; 
		elementArrayZ[i] = 0.0; 
	}

	forceFile.open("out_particleForce.dat",ios::out | ios::binary);
	//aveForceFile.open("aveForce.dat",ios::out | ios::binary);
	//particleDataFile.open("particlePosition.dat",ios::out | ios::binary);
	torqueFile.open("out_torque.txt",ios::out);
	errorFile.open("out_errorlog.txt",ios::out);
	//dataFile.open("out_out.dat",ios::out);
};

Mill::~Mill()
{
	delete id;

	//outputFile.close();
	//dataFile.close();
	forceFile.close();
	torqueFile.close();
	errorFile.close();
	cout<<" Mill Deleted "<<endl;
};

void Mill::run()
{
	time_t start;
	time(&start);
	time_t currTime;
	time(&currTime);
	ofstream outputFile("out_dump",ios::out);
	outputFile<<"START TIME  "<<asctime(localtime(&currTime))<<endl;
	outputFile.close();

	id = new InputData();
	id->readData();
	gapLimit = gapAllowed*id->lengthFactor;
	diskY = 0.0;
	diskThick = id->diskThick;
	diskDia = id->diskDia; 
	shaftDia = id->shaftDia;
	//gForce = id->pMass; //gravitational force
	collisionGapCoeff = collGapCoeff*id->lengthFactor;

   	dx = (id->chamberInnerDia)/(iIndex*1.0);
	dy = (id->chamberInnerDia)/(jIndex*1.0);
	dz = (id->chamberLength)/(kIndex*1.0);  
	cellVolume = dx*dy*dz;
	cout<<" MASS "<<id->pMass<<endl;
	cout<<"dx  "<<dx/id->lengthFactor<<"  dy  "<<dy/id->lengthFactor<<"  dz"<<dz/id->lengthFactor<<endl;
	
	devideIntoCell();
	cout<<" devideIntoCell()"<<endl;
    allocateNeighbouringCell();
    cout<<"allocateNeighbouringCell()"<<endl;
	generateMill();
    cout<<"generateMill()"<<endl;
	startMill();
    cout<<" startMill()"<<endl;
	//updateCfdCells(noOfCfdFrames);//initial condition for CFD cells
	simulate();
    time_t end;
    time(&end);
    double diff = difftime(end, start);
    cout<<"Exe Time  "<<diff<<" s."<<endl;

	cout<<" PARTICLE INSERTION COMPLETE "<<endl;
	time(&currTime);
	ofstream outputFile2("out_dump",ios::app);
	outputFile2<<"END TIME  "<<asctime(localtime(&currTime))<<endl;
	outputFile2.close();
	cout<<"DONE"<<endl;
	system("PAUSE");
};

void Mill::generateMill()
{
	//---------  Disk code ------------/
	generateDisc();
	//---------------------------------

	//---------- Shaft --------------/
	componentType[compIndex] = 2; //type is shaft
	componentR[compIndex] = 0.0; 
	componentDia[compIndex] = id->shaftDia;
	componentTheta[compIndex] = 0.0;
	componentZ[compIndex] = id->chamberLength*0.5;  //half length of the shaft
	componentLength[compIndex] = id->chamberLength;
	componentAngVelocity[compIndex] = -2.0*M_PI*id->rotSpeed/60.0; //mill rotation speed
	compIndex++;
	//------- End of shaft -----------/

	//------- Drum -------------------/
	componentType[compIndex] = 3; //type is drum
	componentR[compIndex] = 0.0; 
	componentDia[compIndex] = id->chamberInnerDia;
	componentTheta[compIndex] = 0.0;
	componentZ[compIndex] = id->chamberLength*0.5;  //half length of the drum
	componentLength[compIndex] = id->chamberLength;
	componentAngVelocity[compIndex] =  0.0; //drum is stationary
	//componentAngVelocity[compIndex] =   -2.0*M_PI*id->rotSpeed/60.0; 
	compIndex++;
	//-------- End of Drum ------------/

	//------ Without boundary conditions -------------------------------
	//-------- First mill wall ----------------
	componentType[compIndex] = 4; //type is drum end wall
	componentR[compIndex] = 0.0; 
	componentDia[compIndex] = id->chamberOuterDia;
	componentTheta[compIndex] = 0.0;
	componentZ[compIndex] = -5.0;
	componentLength[compIndex] = 10.0;
	componentAngVelocity[compIndex] = 0.0; //drum is stationary
    //componentAngVelocity[compIndex] =   -2.0*M_PI*id->rotSpeed/60.0;
	compIndex++;
	//-------- End of left wall ----------------

	//-------- Second mill wall ---------------
	componentType[compIndex] = 5; //type is drum end wall
	componentR[compIndex] = 0.0; 
	componentDia[compIndex] = id->chamberOuterDia;
	componentTheta[compIndex] = 0.0;
	componentZ[compIndex] = id->chamberLength + 5.0;
	componentLength[compIndex] = 10.0;
	componentAngVelocity[compIndex] = 0.0; //drum is stationary
    //componentAngVelocity[compIndex] =   -2.0*M_PI*id->rotSpeed/60.0;
    compIndex++;
	//---------- End of right wall --------------
	//-----------------------------------------------------------------------

    
	cellDx = dx;
	cellDy = dy;
	cellDz = dz;
	
    counter = 0;
	
	fileCounter = 0;
	//--- Wear data ------------------------
    //diskFile.open("in_dsk180_5_6H.dat",ios::in);
    //diskFile.open("in_dsk100_3_5H.dat",ios::in);
	//diskFile.open("in_worn5H.dat",ios::in);
	//diskFile.open("in_worn10H.dat",ios::in);
	//diskFile.open("in_worn15H.dat",ios::in);
	//diskFile.open("in_wornPIN.dat",ios::in);
	//--------------------------------------

    ifstream diskFile("in_dsk180_5_6H.dat",ios::in);
	if (!diskFile)
	{
		cout<<"DISK FILE MISSING!!!!!!!!!!!!!"<<endl;
		exit(1);
	}
	while(fileCounter < diskFileArraySize1)
	{
		diskFile>>tempX>>tempY>>tempZ;
		double tempRR = sqrt(pow(tempX,2)+pow(tempY,2));
		double tTheta = 0.0;
		if(tempX < 0 )
			tTheta = atan(tempY/tempX) + M_PI;
		else
			tTheta = atan(tempY/tempX);
		tTheta = tTheta + rotorAngPosition;
		tempX = tempRR*cos(tTheta);
		tempY = tempRR*sin(tTheta);
		pointArrayX[fileCounter] = tempX*1e-3*id->lengthFactor;
		pointArrayY[fileCounter] = tempY*1e-3*id->lengthFactor;
		pointArrayZ[fileCounter] = tempZ*1e-3*id->lengthFactor + 0.5*id->diskApart;
        fileCounter++;
	}
	fileCounter = 0;
	while(fileCounter < diskFileArraySize2)
	{
		diskFile>>tempX>>tempY>>tempZ;
		elementArrayX[fileCounter] = tempX;
		elementArrayY[fileCounter] = tempY;
		elementArrayZ[fileCounter] = tempZ;
		fileCounter++;
		
	}
	diskFile.close();
//----- REading CFD node file ------------------------------
		int node_counter = 0;
		char ch1 = ' ';
		ifstream nodeFile("export.csv",ios::in);
		if (!nodeFile)
		{
			cout<<" CFD node file missing !!!!!!!!!!"<<endl;
			exit(1);
		}
		nodeFile.ignore(1000,']');
		nodeFile.ignore(1000,']');
		nodeFile.ignore(1000,']');
		nodeFile.ignore(1000,']');
		nodeFile.ignore(1000,']');
		//cfxCellX = new double[cfxCellPoints];       
		//cfxCellY = new double[cfxCellPoints];       
		//cfxCellZ = new double[cfxCellPoints];   
		while(nodeFile)
		{
			double cellX = 0.0;
			double cellY = 0.0;
			double cellZ = 0.0;
			nodeFile>>cellX>>ch>>cellY>>ch1>>cellZ;
			cfxCellX[node_counter] = cellX*id->lengthFactor;
			cfxCellY[node_counter] = cellY*id->lengthFactor;
			cfxCellZ[node_counter] = cellZ*id->lengthFactor;
			node_counter++;
		}
		nodeFile.close();
//---------------------------------------------------------

/*    double cX = 0.0;
    double cY = 0.0;
    double cZ = 0.0;
    double cVelX = 0.0;
    double cVelY = 0.0;
    double cVelZ = 0.0;
	double cPressure = 0.0;
    int cellAIndex = 0;
	int miniCellIndex = 0;
    ifstream cfdFile("in_velocity.dat", ios::in);
	
	while (cfdFile)
    {
		cfdFile>>cX>>cY>>cZ>>cPressure>>cVelX>>cVelY>>cVelZ;
		cfdCellVelX[cellAIndex] = cVelX;
		cfdCellVelY[cellAIndex] = cVelY;
		cfdCellVelZ[cellAIndex] = cVelZ;
		cfdCellPressure[cellAIndex] = cPressure;
		cellAIndex++;

		if(cellAIndex > cfd_iIndex*cfd_jIndex*cfd_kIndex-1)
        {
			break;
        }
    }
    cfdFile.close();

	//time_t pt;
	//time(&pt);
	
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
	*/
	//time_t pEnd;
	//time(&pEnd);
	//double timeDiff = difftime(pEnd,pt);	
	//cout<<"PRESSURE GRADIENT TIME  "<<timeDiff<<endl;	
	/*
	for(int i=0; i<cfd_iIndex*cfd_jIndex*cfd_kIndex; i++)
	{
			errorFile<<setw(15)<<cfdCellCellX[i]*1e3/id->lengthFactor
			<<setw(15)<<cfdCellCellY[i]*1e3/id->lengthFactor
			<<setw(15)<<cfdCellCellZ[i]*1e3/id->lengthFactor
			<<setw(15)<<cfdCellPGradX[i]*id->lengthFactor
			<<setw(15)<<cfdCellPGradY[i]*id->lengthFactor
			<<setw(15)<<cfdCellPGradZ[i]*id->lengthFactor
			<<setw(15)<<cfdCellVelX[i]
			<<setw(15)<<cfdCellVelY[i]
			<<setw(15)<<cfdCellVelZ[i]<<endl;
	}
	exit(0);
	*/
//---------------------------------------------------------
}

/*******************************************
 * This method divide the system into cells
 ********************************************/
void Mill::devideIntoCell()
{
	counter = 0;
	arrCounter = 0;
	double cellY = 0.0;
	double cellX = 0.0;
	double cellZ = 0.0;
	for(int k=0; k<kIndex; k++)
	{
		cellZ = k*dz + dz/2.0;
		for(int j=0; j<jIndex; j++)
		{
			cellY = j*dy + dy/2.0-(id->chamberInnerDia/2.0);
			for(int i=0; i<iIndex; i++)
			{
				int p = i;
				int q = j*iIndex;
				int r = k*iIndex*jIndex;
				cellX = i*dx + dx/2.0-(id->chamberInnerDia/2.0);
				cellCellNo[p+q+r] = counter+1;
				cellCellX[p+q+r] = cellX;
				cellCellY[p+q+r] = cellY;
				cellCellZ[p+q+r] = cellZ;
				counter++;
                                /*if (k==0 && j==0)
                                {
                                        errorFile<<setprecision(4)<<cellX/id->lengthFactor<<" ";
                                }*/

			}
		}
		//errorFile<<endl;
		//errorFile<<cellZ/id->lengthFactor<<" ";
	}
	//errorFile<<endl;
	//exit(0);
	//counter = 0;
	//arrCounter = 0;
	
	double cfdCellY = 0.0;
	double cfdCellX = 0.0;
	double cfdCellZ = 0.0;
	cfdDx = id->chamberInnerDia/(cfd_iIndex*1.0);
	cfdDy = id->chamberInnerDia/(cfd_jIndex*1.0);
	cfdDz = id->chamberLength/(cfd_kIndex*1.0);
	
	for(int k=0; k<cfd_kIndex; k++)
	{
		cfdCellZ = k*cfdDz + cfdDz/2.0;
		for(int j=0; j<cfd_jIndex; j++)
		{
			cfdCellY = j*cfdDy + cfdDy/2.0-(id->chamberInnerDia/2.0);
			for(int i=0; i<cfd_iIndex; i++)
			{
				int p = i;
				int q = j*cfd_iIndex;
				int r = k*cfd_iIndex*cfd_jIndex;
				cfdCellX = i*cfdDx + cfdDx/2.0-(id->chamberInnerDia/2.0);
				cfdCellCellX[p+q+r] = cfdCellX;
				cfdCellCellY[p+q+r] = cfdCellY;
				cfdCellCellZ[p+q+r] = cfdCellZ;
			}
		}
	}
};
/*********************************************
 * Allocate neighbouring cells
 *********************************************/
/*
 //--- For boundary conditions -----------------------------------------------------
void Mill::allocateNeighbouringCell()
{
	for(int k=0; k<kIndex; k++)
	{
		for(int j=0; j<jIndex; j++)
		{
			for(int i=0; i<iIndex; i++)
			{
				int index = 0;
				int p = i;
				int q = j*iIndex;
				int r = k*iIndex*jIndex;
				//neighCellList[index] = cellArray[i][j][k]->cellNo;//cell itself is a neighbour//
				neighCellList[index] = cellCellNo[p+q+r];//cell itself is a neighbour//

				index++;
				if(i!=0)
				{
					if(j!=0)
					{
						if(k!=kIndex-1)
						{	
							//neighCellList[index] = cellArray[i-1][j-1][k+1]->cellNo;
							p = i-1;
							q = (j-1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i-1][j-1][0]->cellNo;
							p = i-1;
							q = (j-1)*iIndex;
							r = 0;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}

						if(k!=0)
						{
							//neighCellList[index] = cellArray[i-1][j-1][k-1]->cellNo;
							p = i-1;
							q = (j-1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i-1][j-1][kIndex-1]->cellNo;
							p = i-1;
							q = (j-1)*iIndex;
							r = (kIndex-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						//neighCellList[index] = cellArray[i-1][j-1][k]->cellNo;
						p = i-1;
						q = (j-1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];

						index++;
					}	
					if(j!=jIndex-1)
					{
						if(k!=kIndex-1)
						{
							//neighCellList[index] = cellArray[i-1][j+1][k+1]->cellNo;
							p = i-1;
							q = (j+1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i-1][j+1][0]->cellNo;
							p = i-1;
							q = (j+1)*iIndex;
							r = 0;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
					
						if(k!=0)
						{
							//neighCellList[index] = cellArray[i-1][j+1][k-1]->cellNo;
							p = i-1;
							q = (j+1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i-1][j+1][kIndex-1]->cellNo;
							p = i-1;
							q = (j+1)*iIndex;
							r = (kIndex-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						//neighCellList[index] = cellArray[i-1][j+1][k]->cellNo;
						p = i-1;
						q = (j+1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}

					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i-1][j][k+1]->cellNo;
						p = i-1;
						q = j*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i-1][j][0]->cellNo;
						p = i-1;
						q = j*iIndex;
						r = 0;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
				
					if(k!=0)
					{
						//neighCellList[index] = cellArray[i-1][j][k-1]->cellNo;
						p = i-1;
						q = j*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i-1][j][kIndex-1]->cellNo;
						p = i-1;
						q = j*iIndex;
						r = (kIndex-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					//neighCellList[index] = cellArray[i-1][j][k]->cellNo;
					p = i-1;
					q = j*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				
				}//end of i!=0

				if(i!=iIndex-1)
				{
					if(j!=0)
					{	
						if(k!=kIndex-1)
						{
							//neighCellList[index] = cellArray[i+1][j-1][k+1]->cellNo;
							p = i+1;
							q = (j-1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i+1][j-1][0]->cellNo;
							p = i+1;
							q = (j-1)*iIndex;
							r = 0;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						if(k!=0)
						{
							//neighCellList[index] = cellArray[i+1][j-1][k-1]->cellNo;
							p = i+1;
							q = (j-1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i+1][j-1][kIndex-1]->cellNo;
							p = i+1;
							q = (j-1)*iIndex;
							r = (kIndex-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						//neighCellList[index] = cellArray[i+1][j-1][k]->cellNo;
						p = i+1;
						q = (j-1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
			
					}
					
					if(j!=jIndex-1)
					{
						if(k!=kIndex-1)
						{
							//neighCellList[index] = cellArray[i+1][j+1][k+1]->cellNo;
							p = i+1;
							q = (j+1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i+1][j+1][0]->cellNo;
							p = i+1;
							q = (j+1)*iIndex;
							r = 0;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						
						if(k!=0)
						{
							//neighCellList[index] = cellArray[i+1][j+1][k-1]->cellNo;
							p = i+1;
							q = (j+1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						else
						{
							//neighCellList[index] = cellArray[i+1][j+1][kIndex-1]->cellNo;
							p = i+1;
							q = (j+1)*iIndex;
							r = (kIndex-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						//neighCellList[index] = cellArray[i+1][j+1][k]->cellNo;
						p = i+1;
						q = (j+1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
			
					}

					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i+1][j][k+1]->cellNo;
						p = i+1;
						q = j*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i+1][j][0]->cellNo;
						p = i+1;
						q = j*iIndex;
						r = 0;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					if(k!=0)
					{
						//neighCellList[index] = cellArray[i+1][j][k-1]->cellNo;
						p = i+1;
						q = j*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i+1][j][kIndex-1]->cellNo;
						p = i+1;
						q = j*iIndex;
						r = (kIndex-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					//neighCellList[index] = cellArray[i+1][j][k]->cellNo;
					p = i+1;
					q = j*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				
				}//end of i!=13

				
				if(j!=0)
				{
					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i][j-1][k+1]->cellNo;
						p = i;
						q = (j-1)*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i][j-1][0]->cellNo;
						p = i;
						q = (j-1)*iIndex;
						r = 0;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}

					if(k!=0)
					{
						//neighCellList[index] = cellArray[i][j-1][k-1]->cellNo;
						p = i;
						q = (j-1)*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i][j-1][kIndex-1]->cellNo;
						p = i;
						q = (j-1)*iIndex;
						r = (kIndex-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					//neighCellList[index] = cellArray[i][j-1][k]->cellNo;
					p = i;
					q = (j-1)*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				}

				if(j!=jIndex-1)
				{
					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i][j+1][k+1]->cellNo;
						p = i;
						q = (j+1)*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i][j+1][0]->cellNo;
						p = i;
						q = (j+1)*iIndex;
						r = 0;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
				
					if(k!=0)
					{
						//neighCellList[index] = cellArray[i][j+1][k-1]->cellNo;
						p = i;
						q = (j+1)*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					else
					{
						//neighCellList[index] = cellArray[i][j+1][kIndex-1]->cellNo;
						p = i;
						q = (j+1)*iIndex;
						r = (kIndex-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					}
					//neighCellList[index] = cellArray[i][j+1][k]->cellNo;
					p = i;
					q = (j+1)*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				}

				if(k!=kIndex-1)
				{
					//neighCellList[index] = cellArray[i][j][k+1]->cellNo;
					p = i;
					q = j*iIndex;
					r = (k+1)*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				}
				else
				{
					//neighCellList[index] = cellArray[i][j][0]->cellNo;
					p = i;
					q = j*iIndex;
					r = 0;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				}
				if(k!=0)
				{
					//neighCellList[index] = cellArray[i][j][k-1]->cellNo;
					p = i;
					q = j*iIndex;
					r = (k-1)*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				}
				else
				{
					//neighCellList[index] = cellArray[i][j][kIndex-1]->cellNo;
					p = i;
					q = j*iIndex;
					r = (kIndex-1)*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
				}
				p = i;
				q = j*iIndex;
				r = k*iIndex*jIndex;
				setNeighbourCell((p+q+r),neighCellList,index);
			}//end of i loop
		}//end of j loop
	}//end of k loop
};

//------------------------------------------------------------------------------------
*/

//----- Without boundary conditions --------------------------------------------------- 
void Mill::allocateNeighbouringCell()
{
	for(int k=0; k<kIndex; k++)
	{
		for(int j=0; j<jIndex; j++)
		{
			for(int i=0; i<iIndex; i++)
			{
				int index = 0;
				
				//neighCellList[index] = cellArray[i][j][k]->cellNo;//cell itself is a neighbour//
				int p = i;
				int q = j*iIndex;
				int r = k*iIndex*jIndex;
				neighCellList[index] = cellCellNo[p+q+r];

				index++;
				if(i!=0)
				{
					if(j!=0)
					{
						if(k!=kIndex-1)
						{	
							//neighCellList[index] = cellArray[i-1][j-1][k+1]->cellNo;
							p = (i-1);
							q = (j-1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}

						//neighCellList[index] = cellArray[i-1][j-1][k]->cellNo;
						p = (i-1);
						q = (j-1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;

						if(k!=0)
						{
							//neighCellList[index] = cellArray[i-1][j-1][k-1]->cellNo;
							p = (i-1);
							q = (j-1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
					}	
					if(j!=jIndex-1)
					{
						if(k!=kIndex-1)
						{
							//neighCellList[index] = cellArray[i-1][j+1][k+1]->cellNo;
							p = (i-1);
							q = (j+1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
						
						//neighCellList[index] = cellArray[i-1][j+1][k]->cellNo;
						p = (i-1);
						q = (j+1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
					
						if(k!=0)
						{
							//neighCellList[index] = cellArray[i-1][j+1][k-1]->cellNo;
							p = (i-1);
							q = (j+1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
						}
					}

					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i-1][j][k+1]->cellNo;
						p = (i-1);
						q = j*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}
					//neighCellList[index] = cellArray[i-1][j][k]->cellNo;
					p = (i-1);
					q = j*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
					

					if(k!=0)
					{
						//neighCellList[index] = cellArray[i-1][j][k-1]->cellNo;
						p = (i-1);
						q = j*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}
				}//end of i!=0

				if(i!=iIndex-1)
				{
					if(j!=0)
					{	
						if(k!=kIndex-1)
						{
							//neighCellList[index] = cellArray[i+1][j-1][k+1]->cellNo;
							p = (i+1);
							q = (j-1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
							
						}

						//neighCellList[index] = cellArray[i+1][j-1][k]->cellNo;
						p = (i+1);
						q = (j-1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						

						if(k!=0)
						{
							//neighCellList[index] = cellArray[i+1][j-1][k-1]->cellNo;
							p = (i+1);
							q = (j-1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
							
						}
			
					}
					
					if(j!=jIndex-1)
					{
						if(k!=kIndex-1)
						{
							//neighCellList[index] = cellArray[i+1][j+1][k+1]->cellNo;
							p = (i+1);
							q = (j+1)*iIndex;
							r = (k+1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
							
						}
						
						//neighCellList[index] = cellArray[i+1][j+1][k]->cellNo;
						p = (i+1);
						q = (j+1)*iIndex;
						r = k*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						

						if(k!=0)
						{
							//neighCellList[index] = cellArray[i+1][j+1][k-1]->cellNo;
							p = (i+1);
							q = (j+1)*iIndex;
							r = (k-1)*iIndex*jIndex;
							neighCellList[index] = cellCellNo[p+q+r];
							index++;
							
						}
					}

					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i+1][j][k+1]->cellNo;
						p = (i+1);
						q = j*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}
					//neighCellList[index] = cellArray[i+1][j][k]->cellNo;
					p = (i+1);
					q = j*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
					

					if(k!=0)
					{
						//neighCellList[index] = cellArray[i+1][j][k-1]->cellNo;
						p = (i+1);
						q = j*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}
				}//end of i!=13

				
				if(j!=0)
				{
				
					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i][j-1][k+1]->cellNo;
						p = i;
						q = (j-1)*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}

					//neighCellList[index] = cellArray[i][j-1][k]->cellNo;
					p = i;
					q = (j-1)*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
					

					if(k!=0)
					{
						//neighCellList[index] = cellArray[i][j-1][k-1]->cellNo;
						p = i;
						q = (j-1)*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}
				}

				if(j!=jIndex-1)
				{
					if(k!=kIndex-1)
					{
						//neighCellList[index] = cellArray[i][j+1][k+1]->cellNo;
						p = i;
						q = (j+1)*iIndex;
						r = (k+1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}
					//neighCellList[index] = cellArray[i][j+1][k]->cellNo;
					p = i;
					q = (j+1)*iIndex;
					r = k*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
					

					if(k!=0)
					{
						//neighCellList[index] = cellArray[i][j+1][k-1]->cellNo;
						p = i;
						q = (j+1)*iIndex;
						r = (k-1)*iIndex*jIndex;
						neighCellList[index] = cellCellNo[p+q+r];
						index++;
						
					}
				}

				if(k!=kIndex-1)
				{
					//neighCellList[index] = cellArray[i][j][k+1]->cellNo;
					p = i;
					q = j*iIndex;
					r = (k+1)*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
					
				}
				if(k!=0)
				{
					//neighCellList[index] = cellArray[i][j][k-1]->cellNo;
					p = i;
					q = j*iIndex;
					r = (k-1)*iIndex*jIndex;
					neighCellList[index] = cellCellNo[p+q+r];
					index++;
					
				}
				p = i;
				q = j*iIndex;
				r = k*iIndex*jIndex;
				setNeighbourCell((p+q+r),neighCellList,index);
			}//end of j loop
		}//end of i loop
	}//end of k loop
};
//----------------------------------------------------------------------------------

