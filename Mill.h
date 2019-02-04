/*****************************************
 * The driving class for the mill
 * author: Chandana Jayasundara 
 * SIMPAS UNSW
 * ver: 5.0
 * May 2007
 *****************************************/
#include<cmath>
#include "InputData.h"
#include<string.h>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<ctime>
#include <iomanip>

#define kIndex 52
#define jIndex 56
#define iIndex 56
#define cfd_iIndex 112
#define cfd_jIndex 112
#define cfd_kIndex 104
#define neighCfdPoints 5
#define cfxCellPoints 323000
#define cfd_neighbourRegion 2*1e-3
#define noOfCfdFrames 36
#define neighCellListSize 27
#define neighParListSize 20 //no of particles in a cell

#define definedRegion 3*1e-3//radius of the spherical region defined for region cells =  diameter of largest particle
#define particleArraySize 225000
#define diskFileArraySize1 3074 // number of finite element data points in the disk file
#define diskFileArraySize2 6168 // number of finite element data connectivity in the disk file
//#define diskFileArraySize1 2668 // number of finite element data points in the disk file
//#define diskFileArraySize2 5352 // number of finite element data connectivity in the disk file

#define holeCenterRadius 60*1e-3
//---- wornHole 1 ---------------
//#define diskFileArraySize1 5377 // number of finite element data points in the disk file
//#define diskFileArraySize2 10770 // number of finite element data connectivity in the disk file
//-------------------------------
/*
//---- wornHole 2 ---------------/
#define diskFileArraySize1 5515 // number of finite element data points in the disk file
#define diskFileArraySize2  11046// number of finite element data connectivity in the disk file
//-------------------------------

//---- wornHole 3 ---------------
#define diskFileArraySize1 5691 // number of finite element data points in the disk file
#define diskFileArraySize2  11398// number of finite element data connectivity in the disk file
//-------------------------------

//---- pin ---------------------
#define diskFileArraySize1 3595 // number of finite element data points in the disk file
#define diskFileArraySize2 7186 // number of finite element data connectivity in the disk file
//-------------------------------
*/

#define lengthFactor1 3.6*1e-3 //by decreasing this factor x-coordinate of region cell can be decreased (outer)  
#define lengthFactor2 6.0*1e-3 //by decreasing this factor x-coordinates of region cell can be increased (inner)
#define cellParListSize 25 //maximum number of particles within a cell 
//#define regionCellParSize 40 //maximum number of particles within a region cell

#define collGapCoeff 0.06*1e-3 //allowed gap for a collision 
//#define maxDisplacement 1.0*1e-3

#define rIn 3*1e-3 //par diameter - largest particle
#define rOut 4.65*1e-3 // rIn*1.55 - largest particle
#define allowedDistance 0.825*1e-3 // (rOut-rIn)/2 - smallest particle

#define compArraySize 22 //no of mill components
#define hollowCompArraySize 6 //number of sub components
//--- Material constants for wear calculation - SAE 1055 steel(fully hardened)----
#define elasticLoadLimit  34.345e9 //g/cm^2;
#define defWearFactor 1.4032e11 //energy required to remove unit volume of material gcm/cm^3
#define KConstant 66.8 //cm/s
#define hardness 1128*1e6 //(Pascal)
#define expKConstant 0.0229 //non dimensional constant

//---- Wear data --------------------------------------------
#define holeOffset 0.0007725 //m
#define holeOffset2 0.002465 //m
#define holeOffsetAngle 2.0 //degrees
#define holeOffsetAngle2 6.0 //degrees

#define zMinimum 0.01 //minimum limit close to disc (m)
#define zMaximum 0.029 //maximmum limit close to disc (m)
#define wearGapLimit 0.004 //(m)
#define radiusAdjustment 0.001 //m
#define noOfHoleDivisions 30 //hole is divided into divisions 
#define alpha 1.2353 //hole r1 and r2
#define alpha2 1.2243 // hole r2 and r3
#define alpha3 1.1812 // hole r3 and r4
//------------------------------------------------------------
#define pointArraySize 20 

#define arraySz 20 //particle contact list
#define collisionIndex 40 //instantanous contact particles
#define neighListArraySz 50 //neighbour list of a particle
#define gapAllowed 0.45*1e-3 //maximum allowed overlap
//#define upperLimit -120
//---- Constants for fluid force --------
//#define dragCoeff 0.1
#define fluidDensity 1300
#define fluidViscosity 0.002


using namespace std;

#ifndef MILL_H
#define MILL_H

class Mill
{
	public:
		Mill();
		~Mill();
		void initiate();
		void devideIntoCell();
		void allocateNeighbouringCell();
		void startMill();
		void run();
		void generateMill();
		bool insertParticle();
		//bool insertParticle(double x,double y,double z);
		bool allocateIntoCells(double pX, double pY, double pZ, int pAIndex);
		void allocateNeighbourList();
		void updateNeighList();
		void simulate();
		void moveParticle();
		void startCalculation();
		void findContact();
		void findContactForce(int contact);
		void calculateForce();
		void generateDisc();
		int vecAdd(double *v1,double *v2);
		int vecSubstract(double *v1,double *v2);
		int unitVector(double *v);
		int sclVecMult(double scl, double *v);
		int projVec(double *v,double *n,int type);
		int crossProduct(double *v1,double *v2);
		int crossProd(double *v1,double *v2);
		double dotProduct(double *v1,double *v2);
		double vecMag(double *v);
		void calculateParForce();
		void calForceBetweenPart();
		void writeDumpFile(double sT, double tT,int cfdFC);
		void writeCfdDetails();
		void writeParticlePosition();
		void writeDiskPosition();
		void writeZonePosition();
		void writePositionData();
        void checkBoundary(int prevB);
		void checkAllBoundary();
		void checkHubBoundary(int prevB);
		void checkDrumBoundary(int prevB);
		void checkShaftBoundary(int prevB);
		void checkSpokeBoundary(int prevB);

		void writeCollisionFile();
		int findSubCompContact(int compI,int subCompI);
		int findCompContact(int compI);
		int findDiskContact(int compI);
		int findShaftContact(int compI);
		int findDrumContact(int compI);
		int findFirstEndWallContact(int compI);
		int findSecondEndWallContact(int compI);
		int findSpokeContact(int compI);
		int findHubContact(int compI);

		double findParVolume(double pD, double defR, double d);

		bool checkData(int showKE,int parI, int arrI, double x,double y,double z,double cF, 
							double iE,double sE, double gap, double kE,double tkE, double tW);
		bool deleteContact(int parI, int arrIndex);
		bool deleteContactNPAR(int parI, int arrIndex);
		void deleteContactList(int parI);
		void setContactList(int pI, int nParI);
		void insertNeighParticle(int parI, int pI);
		void deleteNeighParticle(int parI, int pI);
		void setMiniCellArray(int cellI, double dx, double dy, double dz);
		void setNeighbourCell(int cI, int nL[], int size);
		void deleteParticle(int cI, int pN);
		bool insertParticleToCell(int cI, int pN);

	private:
		char ch;
		int stepNumber,noOfCycles;//number of steps counter used for averagaing CFD data, number of disc rotations
		int frameStepCounter;
		bool firstStep;
		double gapLimit;
		double maxTime,startingTime,totalTime;
		double collisionGapCoeff;
		bool firstInsertion,firstTime,insertable;
		int zoneCounter,counter,tempColor,count,fileCounter,arrCounter,timeCounter;
		int cellParIndex,cellNo;
		double tempX,tempY,tempZ,tempR,dxDot,dyDot,dzDot,xAngVel,yAngVel,zAngVel;
		double dist;
		double rotorSpeed, tempTheta;
		int noOfParticles,index;
		double rotorAngPosition;
		double tempAngPosition;
		double diskAngPosition;
		//time_t currTime;
		int neighCellList[27]; // array for neighbouring cells
		int parList[cellParListSize]; // array for particle list in the cell
		InputData *id;
		
		//--------- Force details ------------------------------
		double instImpactEnergy;
		double nrmForce,nrmDsp,tngDampForce,instRelativeVel;
		double dsmax,dd,fdt,dti,fdtx,fdty,fdtz;
		double instShearForce,instCompForce,instShearEnergy,collEnergy,tangCollEnergy,collPointX,collPointY,collPointZ; 
		double ipX,ipY,ipZ,xx,yy,zz,xy;
		double jpX,jpY,jpZ,nrmCntForce,nrmDampForce;
		double dsmaxCff,tempDist,gForce,torque,totalTorque,totTime;
		double diskThick ,diskDia,shaftDia,pXY;
		double sfCoff, normDampC;
		double r1,r2,r3,ipDia,jpDia,rStar,nrmVel,tangVel;
		double val1,val2,val3,temp,tempVal;
		double dragFX, dragFY, dragFZ, bouyancyF;
		int contact,timeCount,timeCount1,boundary;
        //double totalWear;
        double tempVec[3];
		double vec[3];
		double tipContPointDisp[3];
		double disp[3];
		double contPointDisp[3];
		double unitVec[3];
		double ipRVec[3];
		double jpRVec[3];
		double ijVec[3];
		double rotVel[3];
		double ipCntPntVel[3];
		double jpCntPntVel[3];
		double cntPntVel[3];
		double totForce[3];
		double rotorForce[3];
		double rotorTotForce[3];
		double pMom[3];
		double torqueOnMill[3];
		double momentum[3];
		double momentumR[3];
		double pForce[3];
		double drumRVec[3];
		double shaftRVec[3];
		double diskRVec[3];
		double holeRVec[3];
		double shaftCntPntVel[3];
		double diskCntPntVel[3];
		double shaftAngVel[3];

		//------- Region Cell details ----
		double ddR,ddTheta,ddZ;

		//-------- Cell details -----------
		//double *tempXArray,*tempYArray,*tempZArray;
		double cellVolume;
		int	  cellCellNo[iIndex*jIndex*kIndex];
		double cellCellX[iIndex*jIndex*kIndex];
		double cellCellY[iIndex*jIndex*kIndex];
		double cellCellZ[iIndex*jIndex*kIndex];
		//double cellVelX[iIndex*jIndex*kIndex];
		//double cellVelY[iIndex*jIndex*kIndex];
		//double cellVelZ[iIndex*jIndex*kIndex];
		int  cellNoOfContacts[iIndex*jIndex*kIndex];
		double cellImpactEnergy[iIndex*jIndex*kIndex];
		double cellShearEnergy[iIndex*jIndex*kIndex];
		double cellCollisionEnergy[iIndex*jIndex*kIndex];
		double cellTangCollisionEnergy[iIndex*jIndex*kIndex];

		double cellPorosity[iIndex*jIndex*kIndex];
		//double cellPressure[iIndex*jIndex*kIndex];
		//double cellPGradX[iIndex*jIndex*kIndex];
		//double cellPGradY[iIndex*jIndex*kIndex];
		//double cellPGradZ[iIndex*jIndex*kIndex];

		int cellCellList[iIndex*jIndex*kIndex][neighCellListSize];
		int cellParticleList[iIndex*jIndex*kIndex][neighParListSize];
		int cellCellListSize[iIndex*jIndex*kIndex];
		int cellParticleListSize[iIndex*jIndex*kIndex];

		double cfdCellCellX[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellCellY[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellCellZ[cfd_iIndex*cfd_jIndex*cfd_kIndex];

		double cfdCellVelX[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellVelY[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellVelZ[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellPGradX[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellPGradY[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellPGradZ[cfd_iIndex*cfd_jIndex*cfd_kIndex];
		double cfdCellPressure[cfd_iIndex*cfd_jIndex*cfd_kIndex];

                double pre_cfdCellVelX[cfd_iIndex*cfd_jIndex*cfd_kIndex];
                double pre_cfdCellVelY[cfd_iIndex*cfd_jIndex*cfd_kIndex];
                double pre_cfdCellVelZ[cfd_iIndex*cfd_jIndex*cfd_kIndex];
                double pre_cfdCellPGradX[cfd_iIndex*cfd_jIndex*cfd_kIndex];
                double pre_cfdCellPGradY[cfd_iIndex*cfd_jIndex*cfd_kIndex];
                double pre_cfdCellPGradZ[cfd_iIndex*cfd_jIndex*cfd_kIndex];

		//-- new code -------
		double cfdCellNeighbourDistance[cfd_iIndex*cfd_jIndex*cfd_kIndex][neighCfdPoints];
		double cfdCellNeighbourPressure[cfd_iIndex*cfd_jIndex*cfd_kIndex][neighCfdPoints];
		double cfdCellNeighbourVelX[cfd_iIndex*cfd_jIndex*cfd_kIndex][neighCfdPoints];
		double cfdCellNeighbourVelY[cfd_iIndex*cfd_jIndex*cfd_kIndex][neighCfdPoints];
		double cfdCellNeighbourVelZ[cfd_iIndex*cfd_jIndex*cfd_kIndex][neighCfdPoints];
		int cfdCellNoOfPoints[cfd_iIndex*cfd_jIndex*cfd_kIndex];

		double tempCfxCellX[cfxCellPoints];
		double tempCfxCellY[cfxCellPoints];
		double cfxCellX[cfxCellPoints];
		double cfxCellY[cfxCellPoints];
		double cfxCellZ[cfxCellPoints];
		double cfxCellVelX[cfxCellPoints];
		double cfxCellVelY[cfxCellPoints];
		double cfxCellVelZ[cfxCellPoints];
		double cfxCellPressure[cfxCellPoints];

		int cfdFrameCounter,stepsBetweenFrames;
		double cfdFrameAngle;
		void updateCfdCells(int fcounter);	
		void increaseCollNo();
		//---------- Particle details ---------------
		int parIndex,nParIndex,num;
		double parX,parY,parZ,parDia,gap,collisionGap,nParZ;
		int particleNo[particleArraySize];
		int particleArrayIndex[particleArraySize];
		double particleDiameter[particleArraySize];
		double particleDensity[particleArraySize];
		int particleColor[particleArraySize];
		double particleX[particleArraySize];
		double particleY[particleArraySize];
		double particleZ[particleArraySize];

		double particleVelX[particleArraySize];
		double particleVelY[particleArraySize];
		double particleVelZ[particleArraySize];

		double particleAngVelX[particleArraySize];
		double particleAngVelY[particleArraySize];
		double particleAngVelZ[particleArraySize];
		double particleHisDispX[particleArraySize];
		double particleHisDispY[particleArraySize];
		double particleHisDispZ[particleArraySize];
		double particleDisplacement[particleArraySize];//keeps a record of particle displacement for each time step
		//double particleAverageCompForce[particleArraySize]; 
		//double particleAverageShearForce[particleArraySize];
		//double particleAverageCollE[particleArraySize];
        //double particleInstantCompForce[particleArraySize]; 
		//double particleInstantShearForce[particleArraySize];
       	double particleMass[particleArraySize];
		double particleInertia[particleArraySize];
		bool particleInsertable[particleArraySize];
		bool particleContactCalculated[particleArraySize];
		int particleCellNo[particleArraySize]; //particle belongs to this cell
		int particleNoOfContacts[particleArraySize];//number of contacts for a given time step

		double particleForceX[particleArraySize];
		double particleForceY[particleArraySize];
		double particleForceZ[particleArraySize];
		double particlePForceX[particleArraySize];
		double particlePForceY[particleArraySize];
		double particlePForceZ[particleArraySize];
		double particleDragForceX[particleArraySize];
	    double particleDragForceY[particleArraySize];
        double particleDragForceZ[particleArraySize];
        double particleReynolds[particleArraySize];

		double particleMomentumX[particleArraySize];
		double particleMomentumY[particleArraySize];
		double particleMomentumZ[particleArraySize];
		int totalCollisions[particleArraySize];
		int instantContactList[particleArraySize][collisionIndex];
		int instantContactListSize[particleArraySize];
		int noOfContacts[particleArraySize];
		int neighbourList[particleArraySize][neighListArraySz];
		int neighListSize[particleArraySize];
		int contactList[particleArraySize][arraySz];
		double collisionEnergy[particleArraySize][arraySz];
		double tangCollisoinEnergy[particleArraySize][arraySz];
		double impactEnergyDissipation[particleArraySize][arraySz];
		double shearEnergyDissipation[particleArraySize][arraySz];
		double previousOverlap[particleArraySize][arraySz];
		
		double maxKineticEnergy[particleArraySize][arraySz];
		bool DRUM_OUTER_FACE[particleArraySize];
		bool SHAFT[particleArraySize];
		bool DRUM_SIDE_FACE[particleArraySize];
		bool DISK_SIDE_FACE[particleArraySize];
		bool HUB_OUTER_FACE[particleArraySize];
		bool HUB_SIDE_FACE[particleArraySize];

		double totalWear[particleArraySize];
		double holeWearData1[noOfHoleDivisions];
		double holeWearData2[noOfHoleDivisions];
		double totalWear1,totalWear2;
		//------------------------------------------------------

		//-------- Component details --------------------------
		double 	componentX[compArraySize];
		double 	componentY[compArraySize]; 
		double 	componentZ[compArraySize];
		double  componentR[compArraySize];
		double  componentDia[compArraySize];
		double  componentTheta[compArraySize];
		double 	componentLength[compArraySize];
		double	componentHeight[compArraySize];
		double	componentWidth[compArraySize]; 

		double componentAngVelocity[compArraySize];
		double componentSlidingFriction[compArraySize];
		int noOfSubComponents[compArraySize];
		int componentType[compArraySize]; //1=disk,2=shaft,3=drum,4=drum side wall
		double hollowCompR[compArraySize][hollowCompArraySize];
		double hollowCompTheta[compArraySize][hollowCompArraySize];
		double hollowCompZ[compArraySize][hollowCompArraySize];
		double hollowCompLength[compArraySize][hollowCompArraySize];
		double hollowCompDia[compArraySize][hollowCompArraySize];
		double hollowCompAngVelocity[compArraySize][hollowCompArraySize];
		int compIndex; //array index for number of components

		bool isInside,calculateTorque;
		int isContacted, nearestSpokeIndex;
		double parXDash,parYDash,parZDash;
		double nearestHoleX,nearestHoleY,nearestDiskZ,nearestHubZ; //use in checkBoundary() function
		int shaftIndex, drumIndex, firstEndWallIndex, secondEndWallIndex;
		int compType;
		double hubDia, hubThick;
		double torqueRadius; //Distance to the contact point when calculating torque

		//----------- Arrays used for reading "disk file" -------------
		double pointArrayX[diskFileArraySize1];
		double pointArrayY[diskFileArraySize1];
		double pointArrayZ[diskFileArraySize1];
		double pointData1[diskFileArraySize1];
		double pointData2[diskFileArraySize1];
		double elementArrayX[diskFileArraySize2];
		double elementArrayY[diskFileArraySize2];
		double elementArrayZ[diskFileArraySize2];
		//-------------------------------------------------------------
		
		double xCoordinate,yCoordinate,zCoordinate,pX,pY,pZ,dx,dy,dz,cfdDx,cfdDy,cfdDz,cellDx,cellDy,cellDz,diskY;
		ofstream torqueFile,errorFile,dataFile;
		fstream forceFile;
};

#endif
