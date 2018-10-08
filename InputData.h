#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<cmath>
using namespace std;
#ifndef  INPUTDATA_H
#define INPUTDATA_H
#define dataArraySize 35
#define largestParDia 3.0
#define largestParDensity 2500.0
#define gravitationalAcc 9.81

class InputData
{
	public:
	InputData(); 
	~InputData();
	void readData();
	int noOfParticles,setCycles;;
	int noOfDisks;
	double diskDia;
	double diskThick;
	double diskApart;
	double shaftDia;
	double shaftLength;
	double chamberInnerDia;
	double chamberOuterDia;
	double chamberLength;
	double holeDia;
	int noOfHoles;
	double particleDia,rotSpeed;
	double timeStep,maxTime,maxRotAngle;
	double pElasticMod,ppNormDampC,pwNormDampC,pdNormDampC,pDensity,sfricp,sfricpWall,sfricpDisk,
		rollfric,maxGap,pMass,pInertia,pPoisonR,pYoungMod;
	double dataValue[dataArraySize];
	double refLength, refDensity;
	double lengthFactor,volumeFactor, massFactor,timeFactor,densityFactor,forceFactor,pressureFactor,StressFactor,
		energyFactor,momentFactor,powerFactor,velocityFactor,accFactor,angVelFactor,angAccFactor,freqFactor,inertiaFactor;
	
	
};
#endif		
