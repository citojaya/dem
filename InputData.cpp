
#define PI 3.14159
#include "InputData.h"
InputData::InputData()
{
	setCycles = 0;
	noOfParticles = 0;
	noOfDisks = 0;
	diskDia = 0.0;
	diskThick = 0.0;
	diskApart = 0.0;
	shaftDia = 0.0;
	shaftLength = 0.0;
	chamberInnerDia = 0.0;
	chamberOuterDia = 0.0;
	chamberLength = 0.0;
	holeDia = 0.0;
	noOfHoles = 0;
	particleDia = 0.0;
	rotSpeed = 0.0;
	rollfric = 0.0;
	timeStep = 0.0;
	maxTime = 0.0;
	pElasticMod = 0.0;
	pYoungMod = 0.0;
	ppNormDampC = 0.0;
	pwNormDampC = 0.0;
	pdNormDampC = 0.0;
	pDensity = 0.0;
	sfricp = 0.0;
	sfricpWall = 0.0;
	sfricpDisk = 0.0;
	maxGap = 0.0;
	pMass = 0.0;
	pInertia = 0.0;
	pPoisonR = 0.0;
	maxRotAngle = 0.0;
	for (int i=0;i<dataArraySize ;i++ )
	{
		dataValue[i] = 0.0;
	}
	lengthFactor = 0.0;
	volumeFactor = 0.0;
	massFactor = 0.0;
	timeFactor = 0.0;
	densityFactor = 0.0;
	forceFactor = 0.0;
	pressureFactor = 0.0;
	StressFactor = 0.0;
	energyFactor = 0.0;
	momentFactor = 0.0;
	powerFactor = 0.0;
	velocityFactor = 0.0;
	accFactor = 0.0;
	angVelFactor = 0.0;
	angAccFactor = 0.0;
	freqFactor = 0.0;
};

InputData::~InputData()
{
	//cout<<"END INPUT_DATA"<<endl;
	//delete this;
};

void InputData::readData()
{
	cout<<"START MILL\n";
	
	char name[10];
	int count = 0;
	ifstream File;
	File.open("in_input_data.txt",ios::in);
	
	if(!File)
	{
		cerr<<"Error reading file";
	}
	while(!File.eof())
	{
		File.ignore(10000,'<');
		File.getline(name,10,'>');
		dataValue[count] = atof(name);
		count++;
	}
	File.close();
	
	chamberOuterDia = dataValue[0]*1e-3;	
	chamberInnerDia = dataValue[1]*1e-3;	
	diskDia = dataValue[2]*1e-3;	
	diskThick = dataValue[3]*1e-3;	
	diskApart = dataValue[4]*1e-3;	
	noOfDisks = (int)dataValue[5];	
	noOfHoles = (int)dataValue[6];	
	holeDia = dataValue[7]*1e-3;
	shaftDia = dataValue[8]*1e-3;
	noOfParticles = (int)dataValue[9];	
	particleDia = dataValue[10]*1e-3;
	pDensity = dataValue[11];
	maxGap = dataValue[14]*1e-3;
	
	pYoungMod = dataValue[15];
	pPoisonR = dataValue[17];
	ppNormDampC = dataValue[19];
	pwNormDampC = dataValue[20];
	pdNormDampC = dataValue[21];
	sfricp = dataValue[25];
	sfricpWall = dataValue[26];
	sfricpDisk = dataValue[27];
	rollfric = dataValue[28];
	rotSpeed = dataValue[29];
	timeStep = dataValue[30];
	maxTime = dataValue[31];
	setCycles = (int)dataValue[32];
	pElasticMod = pYoungMod/(1-pow(pPoisonR,2));

	cout<<" E MOD "<<pElasticMod;
	//--- With boundary conditions -----------
	//chamberLength = (diskThick + diskApart)*noOfDisks;	
	//----------------------------------------

	//--- Without boundary conditions -------------------------------
	//chamberLength = noOfDisks*(diskThick + diskApart) + diskApart;	
	chamberLength = 3.0*(diskThick + diskApart);	
	//---------------------------------------------------------------

	cout<<" CHAMBER LENGTH "<<chamberLength;
	pMass = (4.0/3.0)*M_PI*pow((0.5*particleDia),3.0)*pDensity;
	cout<<" PMASS "<<pMass;
	pInertia = 2*pMass*pow(0.5*particleDia,2)/5.0;//in kgm^2//
	cout<<" INERTIA "<<pInertia;

	//---- Scale factors for reduced units ------------
	refLength = largestParDia*1e-3;
	refDensity = largestParDensity;
	lengthFactor = 1.0/refLength;
	volumeFactor = pow(lengthFactor,3);
	massFactor = 6.0/(M_PI*pow(refLength,3)*refDensity);
	timeFactor = sqrt(gravitationalAcc/refLength);
	densityFactor = 6.0/(M_PI*refDensity);
	forceFactor = 6.0/(gravitationalAcc*M_PI*pow(refLength,3)*refDensity);
	pressureFactor = 6.0/(gravitationalAcc*M_PI*refLength*refDensity);
	StressFactor = pressureFactor;
	energyFactor = 6.0/(gravitationalAcc*M_PI*pow(refLength,4)*refDensity);
	momentFactor = energyFactor;
	powerFactor = 6.0/(pow(gravitationalAcc,1.5)*M_PI*pow((double)refLength,3.5)*refDensity);
	velocityFactor = 1.0/sqrt(refLength*gravitationalAcc);
	accFactor = 1.0/gravitationalAcc;
	angVelFactor = sqrt(refLength/gravitationalAcc);
	angAccFactor = refLength/gravitationalAcc;
	freqFactor = sqrt(refLength/gravitationalAcc);
	inertiaFactor = 6.0/(M_PI*pow(refLength,5)*refDensity);
	//-------------------------------------------------
	//------- Reduced units ---------------------------
	chamberOuterDia = chamberOuterDia*lengthFactor;	
	chamberInnerDia = chamberInnerDia*lengthFactor;	
	diskDia			= diskDia*lengthFactor;	
	diskThick		= diskThick*lengthFactor;	
	diskApart		= diskApart*lengthFactor;	
	holeDia			= holeDia*lengthFactor;
	shaftDia		= shaftDia*lengthFactor;
	particleDia		= particleDia*lengthFactor;
	pDensity		= pDensity*densityFactor;
	maxGap			= maxGap*lengthFactor;
					                  
	pYoungMod		= pYoungMod*pressureFactor;
	rotSpeed		= rotSpeed*freqFactor;
	timeStep		= timeStep*timeFactor;
	maxTime			= maxTime*timeFactor;
	pElasticMod		= pYoungMod/(1-pow(pPoisonR,2));

	chamberLength	= chamberLength*lengthFactor;	
	pMass			= pMass*massFactor;
	pInertia		= pInertia*inertiaFactor;
	//-------------------------------------------------
	cout<<" MAss  "<<pMass<<"  densFactor "<<densityFactor<<endl;
	cout<<" Vel Factor "<<velocityFactor<<endl;
};
