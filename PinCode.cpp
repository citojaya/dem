#include "Mill.h"

void Mill::generateDisc()
{
	int hollowCompIndex = 0;
	int noOfSpokes = 4; //each disk has 4 spokes

	//--------  Disks ------------
	for (int j=0; j<noOfSpokes ;j++ )
	{
		noOfSubComponents[compIndex] = hollowCompIndex;
		componentType[compIndex] = 1; //type is rectangle
		//componentR[compIndex] = 32.5; 
		componentR[compIndex] = 30.0*1e-3*id->lengthFactor; 
		componentDia[compIndex] = 0.0;
		componentTheta[compIndex] = j*90.0*M_PI/180.0;
		componentZ[compIndex] = id->diskApart*0.5+id->diskThick*0.5;
		componentLength[compIndex] = 9.0*1e-3*id->lengthFactor;
		componentHeight[compIndex] = 25.0*1e-3*id->lengthFactor;
		componentWidth[compIndex] = 40.0*1e-3*id->lengthFactor;
		componentAngVelocity[compIndex] = 0.0; //mill rotation speed
		compIndex++;
	}

	//----------- End of disks -----------/

	//-------- Spoke hub ----------/
	for (int i=0; i<id->noOfDisks; i++ )
	{
		componentType[compIndex] = 6; //type is hub
		componentR[compIndex] = 0.0; 
		componentDia[compIndex] = 40.0*1e-3*id->lengthFactor;
		hubDia  = 40.0*1e-3*id->lengthFactor;
		componentTheta[compIndex] = 0.0;
		componentZ[compIndex] = id->diskApart*0.5+id->diskThick*0.5;
		componentLength[compIndex] = 9.0*1e-3*id->lengthFactor;
		hubThick = 9.0*id->lengthFactor;
		//componentHeight[compIndex] = 25.0;
		//componentWidth[compIndex] = 9.0;
		componentAngVelocity[compIndex] = -2.0*M_PI*id->rotSpeed/60.0; //mill rotation speed
		compIndex++;
	}
};

int Mill::findCompContact(int compI)
{
	calculateTorque = false; //reset torque calculation condition
	int type = componentType[compI];
	compType = 0;
	switch(type)
	{
		case 1:
			findSpokeContact(compI);
		break;

		case 2:
			findShaftContact(compI);
		break;

		case 3:
			findDrumContact(compI);
		break;

		case 4:
			findFirstEndWallContact(compI);
		break;

		case 5:
			findSecondEndWallContact(compI);
		break;

		case 6:
			findHubContact(compI);
		break;

	}
    return 1;
};
 int Mill::findSpokeContact(int compI)
 {
	//cout<<" SPOKE  "<<endl;
	pXY = sqrt(pow(parX,2)+pow(parY,2));
	calculateTorque = true; //torque calculation is required 
 	contact = 0;
    rStar = 0.5*parDia;
	isContacted = 0; //reset
	sfCoff = id->sfricpDisk;
	float ijVecX = 0.0;
	float ijVecY = 0.0;
	float ijVecZ = 0.0;
	float ijXYZ = 0.0;
	
	float node1Z = 0.0;
	float node2Z = componentLength[compI];
	float theta = componentTheta[compI];
	float ro = componentR[compI];
	float compLength = componentLength[compI];

	//--------- new coordinates of the particle in x'y'z' coordinate system
	parXDash = parX*cos(theta)+parY*sin(theta)- ro;
	parYDash = parY*cos(theta)-parX*sin(theta);
	parZDash = parZ - (componentZ[compI]-compLength*0.5);

	//cout<<" PAR X Y Z"<<parXDash<<"  "<<parYDash<<"   "<<parZDash<<endl;
	nearestSpokeIndex = compI;
	if (parYDash > componentHeight[compI]*0.5 && parYDash < (componentHeight[compI] + parDia + 4*collisionGapCoeff)*0.5 )
	{
		if (parXDash <= componentWidth[compI]*0.5 && parXDash >=  - componentWidth[compI]*0.5)
		{ 
			if (parZDash >= 0 && parZDash <= componentLength[compI])
			{
                gap = fabs(parYDash) - 0.5*parDia - componentHeight[compI]*0.5;
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Top plane GAP "<<gap<<endl;
					}
					unitVec[0] = 0.0;                                                                                                                                                                                                  
					unitVec[1] = 1.0;                                                                                                                                                                                                  
					unitVec[2] = 0.0;                                                                                                                                                                                                 
					//unit vectors for the particle                                                                                                                                                                                    
					ipRVec[0] = -0.5*parDia*(*unitVec);                                                                                                                                                                                
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));                                                                                                                                                                            
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));                                                                                                                                                                            
									
					jpRVec[0] = 0.0;                                                                                                                                                                                              
					jpRVec[1] = componentHeight[compI]*0.5;                                                                                                                                                                                              
					jpRVec[2] = 0.0;                                                                                                                                                                                                   
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta) + ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta) + ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					boundary = -1256;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];
					torqueRadius = 0.0; //torque due to shear force not necessary
					
					findContactForce(1); 
					
					if (HUB_OUTER_FACE[parIndex])                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
						checkHubBoundary(-44);                                                                                                                                                                                             
					}                                                                                                                                                                                                                
																																																										
					isContacted = 1;                                                                                                                                                                                                   
					return isContacted;                                                                                                                                                                                     
				}
				checkSpokeBoundary(-1256);
			}
			else if (parZDash < 0 )
			{
				//cout<<" EDGE EF "<<endl;
				contact = 56;
			}
			else if (parZDash > componentLength[compI])
			{
				//cout<<" EDGE AB "<<endl;
				contact = 12;
			}
		}

		else if (parZDash <= componentLength[compI] && parZDash >= 0)
		{
			if (parXDash < - componentWidth[compI]*0.5)
			{
				//cout<<" EDGE AE" <<endl;
			}
			else if (parXDash > componentWidth[compI]*0.5)
			{
				//errorFile<<" EDGE BF TOP PLANE" <<endl;
				contact = 26;
			}
		}

		else if (parXDash < - componentWidth[compI]*0.5)
		{
			if (parZDash > componentLength[compI])
			{
				//cout<<" CORNER A"<<endl;
			}
			else if (parZDash < 0)
			{
				//cout<<" CORNER E"<<endl;
			}
		}
		else if (parXDash > componentWidth[compI]*0.5)
		{
			if (parZDash > componentLength[compI])
			{
				//cout<<" CORNER  B"<<endl;
				contact = 2;
			}
			else if (parZDash < 0)
			{
				//cout<<" CORNER F"<<endl;
				contact = 6;
			}
		}
	}

	if (parYDash < - componentHeight[compI]*0.5 && parYDash > -(componentHeight[compI] + parDia + 4*collisionGapCoeff)*0.5)
	{
		if (parXDash <= componentWidth[compI]*0.5 && parXDash >=  - componentWidth[compI]*0.5)
		{ 
			if (parZDash >= 0 && parZDash <= componentLength[compI])
			{
				//cout<<" BOTTOM PLANE "<<parXDash<<"   "<<parYDash<<"   "<<parZDash<<endl;
				//system("PAUSE");
				
				gap = fabs(parYDash) - 0.5*parDia - componentHeight[compI]*0.5;
				if (gap < 0)
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Bottom plane GAP "<<gap<<endl;
					}
	
					unitVec[0] = 0.0;                                                                                                                                                                                                  
					unitVec[1] = -1.0;                                                                                                                                                                                                  
					unitVec[2] = 0.0;                                                                                                                                                                                                 
					//unit vectors for the particle                                                                                                                                                                                    
					ipRVec[0] = -0.5*parDia*(*unitVec);                                                                                                                                                                                
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));                                                                                                                                                                            
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));                                                                                                                                                                            
									
					jpRVec[0] = 0.0;                                                                                                                                                                                              
					jpRVec[1] = - componentHeight[compI]*0.5;                                                                                                                                                                                              
					jpRVec[2] = 0.0;                                                                                                                                                                                                   
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta) + ro*cos(theta);                                                                                                                                                                        
					jpRVec[1] = r2*cos(theta) + r1*sin(theta) + ro*sin(theta);                                                                                                                                                                        
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					boundary = -3478; 
					torqueRadius = 0.0; //torque due to shear force not necessary
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					
					if (HUB_OUTER_FACE[parIndex])                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkHubBoundary(-44);                                                                                                                                                                                             
					}                                                                                                                                                                                                                
																																																										
					isContacted = 1;
					return isContacted;                                                                                                                                                                                                  
				}
				checkSpokeBoundary(-3478);
				//cout<<" GAP "<< fabs(parYDash) - 0.5*parDia - componentHeight[compI]*0.5<<endl;
				//system("PAUSE");
			}
			else if (parZDash < 0 )
			{
				//cout<<" EDGE HG "<<endl;
				contact = 78;
			}
			else if (parZDash > componentLength[compI])
			{
				//cout<<" EDGE DC "<<endl;
				contact = 34;
			}
		}

		else if (parZDash <= componentLength[compI] && parZDash >= 0)
		{
			if (parXDash < - componentWidth[compI]*0.5)
			{
				//cout<<" EDGE DH" <<endl;
			}
			else if (parXDash > componentWidth[compI]*0.5)
			{
				//cout<<" EDGE GC" <<endl;
				//errorFile<<" edge GC - BOTTOM "<<endl;
				contact = 37;
			}
		}

		else if (parXDash < - componentWidth[compI]*0.5)
		{
			if (parZDash > componentLength[compI])
			{
				//cout<<" CORNER D"<<endl;
			}
			else if (parZDash < 0)
			{
				//cout<<" CORNER H"<<endl;
			}
		}
		else if (parXDash > componentWidth[compI]*0.5)
		{
			if (parZDash > componentLength[compI])
			{
				//cout<<" CORNER  C"<<endl;
				contact = 3;

			}
			else if (parZDash < 0)
			{
				//cout<<" CORNER G"<<endl;
				contact = 7;

			}
		}
	}// End of bottom plane

	if (parZDash > componentLength[compI] && parZDash < componentLength[compI] + parDia*0.5 + 2*collisionGapCoeff)
	{
		if (parXDash <= componentWidth[compI]*0.5 && parXDash >=  - componentWidth[compI]*0.5)
		{ 
			if (parYDash >= -componentHeight[compI]*0.5 && parYDash <= componentHeight[compI]*0.5)
			{
				//cout<<" ABCD PLANE "<<endl;
 				//cout<<" GAP "<< fabs(parYDash) - 0.5*parDia - componentHeight[compI]*0.5<<endl;
				gap = fabs(parZDash)-node2Z - parDia*0.5;
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" ABCD plane GAP "<<gap<<endl;
					}
	
					unitVec[0] = 0.0;                                                                                                                                                                                                  
					unitVec[1] = 0.0;                                                                                                                                                                                                  
					unitVec[2] = 1.0;                                                                                                                                                                                                 
					//unit vectors for the particle                                                                                                                                                                                    
					ipRVec[0] = -0.5*parDia*(*unitVec);                                                                                                                                                                                
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));                                                                                                                                                                            
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));                                                                                                                                                                            
									
					jpRVec[0] = parXDash;                                                                                                                                                                            
					jpRVec[1] = parYDash;                                                                                                                                                                            
					jpRVec[2] = 0.0;                                                                                                                                                                                                        
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                        
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                        
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					boundary = -1234;
					torqueRadius = pXY; //torque due to shear force
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
																																																										
					isContacted = 1;  
					return isContacted;
					//rotatingSpeed = rotSpeed;                                                                                                                                                                                        
				}
				checkSpokeBoundary(-1234);
				//system("PAUSE");
			}

			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//edge DC
				contact = 34;
			}

			else if (parYDash > componentHeight[compI]*0.5)
			{
				//edge AB
				contact = 12;
			}
		}
		else if (parYDash <= componentHeight[compI]*0.5 && parYDash >= -componentHeight[compI]*0.5)
		{
			if (parXDash < - componentWidth[compI]*0.5)
			{
					//edge AD
			}
			else if (parXDash > componentWidth[compI]*0.5)
			{
					//cout<<"  edge BC  "<<endl;
					contact  = 23;
			}
			
		}

		else if (parXDash < - componentWidth[compI]*0.5)
		{
			if (parYDash > componentHeight[compI]*0.5)
			{
				//cout<<" CORNER A"<<endl;
			}
			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//cout<<" CORNER D"<<endl;
			}
		}
		else if (parXDash > componentWidth[compI]*0.5)
		{
			if (parYDash > componentHeight[compI]*0.5)
			{
				//cout<<" CORNER  B"<<endl;
				contact = 2;
			}
			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//cout<<" CORNER C"<<endl;
				contact = 3;
			}
		}
	}//End of ABCD plane

	if (parZDash < 0 && parZDash > -(parDia*0.5 + 2*collisionGapCoeff))
	{
		if (parXDash <= componentWidth[compI]*0.5 && parXDash >=  - componentWidth[compI]*0.5)
		{ 
			if (parYDash >= -componentHeight[compI]*0.5 && parYDash <= componentHeight[compI]*0.5)
			{
				//cout<<" EFGH PLANE "<<endl;
				gap = fabs(parZDash)-node1Z - parDia*0.5;
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" EFGH plane GAP "<<gap<<endl;
					}
	
					//system("PAUSE");
					unitVec[0] = 0.0;                                                                                                                                                                                                  
					unitVec[1] = 0.0;                                                                                                                                                                                                  
					unitVec[2] = -1.0;                                                                                                                                                                                                 
					//unit vectors for the particle                                                                                                                                                                                    
					ipRVec[0] = -0.5*parDia*(*unitVec);                                                                                                                                                                                
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));                                                                                                                                                                            
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));                                                                                                                                                                            
									
					jpRVec[0] = parXDash;                                                                                                                                                                            
					jpRVec[1] = parYDash;                                                                                                                                                                            
					jpRVec[2] = 0.0;                                                                                                                                                                                                        
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                        
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                        
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					boundary = -5678;  
					torqueRadius = pXY; //torque due to shear force 
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
																																																										
					isContacted = 1;   
					return isContacted;
					//rotatingSpeed = rotSpeed;                                                                                                                                                                                        
				}
				checkSpokeBoundary(-5678);
			}
			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//edge GH
				contact = 78;
			}
			else if (parYDash > componentHeight[compI]*0.5)
			{
				//edge EF
				contact = 56;
			}
		
		}
		else if (parYDash <= componentHeight[compI]*0.5 && parYDash >= -componentHeight[compI]*0.5 )
		{
			if (parXDash < - componentWidth[compI]*0.5)
			{
					//edge EH
				
			}
			else if (parXDash > componentWidth[compI]*0.5)
			{
					//cout<< "edge GF"<<endl;
					contact = 67;
			}
		}
	
		else if (parXDash < - componentWidth[compI]*0.5)
		{
			if (parYDash > componentHeight[compI]*0.5)
			{
				//cout<<" CORNER E"<<endl;
			}
			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//cout<<" CORNER H"<<endl;
			}
		}
		else if (parXDash > componentWidth[compI]*0.5)
		{
			if (parYDash > componentHeight[compI]*0.5)
			{
				//cout<<" CORNER  F"<<endl;
				contact = 6;
			}
			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//cout<<" CORNER G"<<endl;
				contact = 7;
			}
		}
	}//End of EFGH plane

	if (parXDash > componentWidth[compI]*0.5 && parXDash < (componentWidth[compI]+parDia+4*collisionGapCoeff)*0.5)
	{
		if (parZDash >= 0 && parZDash <= componentLength[compI])
		{ 
			if (parYDash <= componentHeight[compI]*0.5 && parYDash >= -componentHeight[compI]*0.5)
			{
				//cout<<" BCGF PLANE "<<endl;
				gap = fabs(parXDash) - 0.5*parDia - componentWidth[compI]*0.5;
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" BCGF plane GAP "<<gap<<endl;
					}
					unitVec[0] = 1.0;                                                                                                                                                                                                  
					unitVec[1] = 0.0;                                                                                                                                                                                                  
					unitVec[2] = 0.0;                                                                                                                                                                                                 
					//unit vectors for the particle                                                                                                                                                                                    
					ipRVec[0] = -0.5*parDia*(*unitVec);                                                                                                                                                                                
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));                                                                                                                                                                            
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));                                                                                                                                                                            
									
					jpRVec[0] = componentWidth[compI]*0.5;                                                                                                                                                                                             
					jpRVec[1] = 0.0;                                                                                                                                                                                              
					jpRVec[2] = 0.0;                                                                                                                                                                                                   
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                        
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                        
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					boundary = -2367;   
					torqueRadius = pXY - 0.5*parDia; //torque due to shear force
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
																																																										
					isContacted = 1;  
					return isContacted;
					//rotatingSpeed = rotSpeed;                                                                                                                                                                                        
				}
				checkSpokeBoundary(-2367);
	
 				//system("PAUSE");
			}

			else if (parYDash > componentHeight[compI]*0.5)
			{
				//edge BF
				//errorFile<<"BCGF"<<endl;
				contact = 26;
			}

			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//edge GC
				//errorFile<<" edge GC - BCGF"<<endl;
				contact = 37;
			}
		}
		else if (parYDash <= componentHeight[compI]*0.5 && parYDash >= -componentHeight[compI]*0.5)
		{
			if (parZDash < 0)
			{
				//edge GF
				contact = 67;
			}
			else if (parZDash > node2Z)
			{
				//cout<<"  edge BC  "<<endl;
				contact  = 23;
			}
		}


		else if (parZDash < 0)
		{
			if (parYDash > componentHeight[compI]*0.5)
			{
				//cout<<" CORNER F"<<endl;
				contact = 6;
			}
			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//cout<<" CORNER G"<<endl;
				contact = 7;
			}
		}
		else if (parZDash > componentLength[compI]*0.5)
		{
			if (parYDash > componentHeight[compI]*0.5)
			{
				//cout<<" CORNER  B"<<endl;
				contact = 2;
			}
			else if (parYDash < -componentHeight[compI]*0.5)
			{
				//cout<<" CORNER C"<<endl;
				contact = 3;
			}
		}
	}//End of BCGF plane

	switch(contact)
	{
		case 2: //corner B

				jpX = componentWidth[compI]*0.5;                                                             
				jpY = componentHeight[compI]*0.5;                                                          
				jpZ = node2Z;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					//system("PAUSE");
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Corner B GAP "<<gap<<endl;
					}
					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);
					jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
					jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
					torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force
					boundary = -2;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					return isContacted;
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;                                                                                                                                                                                                   
				}
				checkSpokeBoundary(-2);
		break; 

		case 3: //corner C
				jpX = componentWidth[compI]*0.5;                                                             
				jpY = -componentHeight[compI]*0.5;                                                          
				jpZ = node2Z;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					//system("PAUSE");
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Corner C GAP "<<gap<<endl;
					}

					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);   
					jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
					jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
					torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force
					boundary = -3;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;  
					return isContacted;
				}
				checkSpokeBoundary(-3);
	
		break;

		case 6: //corner F
				jpX = componentWidth[compI]*0.5;                                                             
				jpY = componentHeight[compI]*0.5;                                                          
				jpZ = node1Z;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					//system("PAUSE");
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Corner F GAP "<<gap<<endl;
					}
					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta); 
					jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
					jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
					torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force
					boundary = -6;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;    
					return isContacted;
				}
				checkSpokeBoundary(-6);
		break;

		case 7: //corner G
				jpX = componentWidth[compI]*0.5;                                                             
				jpY = -componentHeight[compI]*0.5;                                                          
				jpZ = node1Z;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Corner G GAP "<<gap<<endl;
					}

					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);     
					jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
					jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
					torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force
					boundary = -7;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;    
					return isContacted;
				}
				checkSpokeBoundary(-7);
		break;

		case 23: //edge BC
					//system("PAUSE");
					jpX = componentWidth[compI]*0.5;                                                             
					jpY = parYDash;                                                          
					jpZ = node2Z;                                                                          
					gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
					//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
					//cout<<" GAP "<<gap<<endl;
					if (gap < 0 )
					{
						if (gap < -gapLimit)
						{
							errorFile<<" SPOKE "<<compI<<" Edge BC GAP "<<gap<<endl;
						}
		
						ijVecX = parXDash - jpX;
						ijVecY = parYDash - jpY;
						ijVecZ = parZDash - jpZ;
						ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

						unitVec[0] = ijVecX/ijXYZ;
						unitVec[1] = ijVecY/ijXYZ;
						unitVec[2] = ijVecZ/ijXYZ;
										//unit vectors for the particle
						ipRVec[0] = -0.5*parDia*(*unitVec);
						ipRVec[1] = -0.5*parDia*(*(unitVec+1));
						ipRVec[2] = -0.5*parDia*(*(unitVec+2));
										//unit vectors for the hole contact point
						jpRVec[0] = jpX;     
						jpRVec[1] = jpY; 
						jpRVec[2] = 0.0; 
																																																														   
						// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
						r1 = ipRVec[0];                                                                                                                                                                                                    
						r2 = ipRVec[1];                                                                                                                                                                                                    
						//r3 = ipRVec[2];                                                                                                                                                                                                  
						ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
						ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
						//ipRVec[2] = does not change                                                                                                                                                                                      
																																																														   
						r1 = jpRVec[0];                                                                                                                                                                                                    
						r2 = jpRVec[1];                                                                                                                                                                                                    
						//r3 = jpRVec[2];                                                                                                                                                                                                  
						jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
						jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
						//jpRVec[2] = does not change                                                                                                                                                                                      
																																																														   
						r1 = unitVec[0];                                                                                                                                                                                                   
						r2 = unitVec[1];                                                                                                                                                                                                   
						//r3 = unitVec[2];                                                                                                                                                                                                 
						unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
						unitVec[1] = r2*cos(theta) + r1*sin(theta);  
						jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
						jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
						torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force
						boundary = -23;                                                                                                                                                                                                     
						shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
						findContactForce(1); 
						/*
						if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
						{                                                                                                                                                                                                                  
								checkBoundary(-2);                                                                                                                                                                                             
						}*/                                                                                                                                                                                                                 
						isContacted = 1; 
						return isContacted;
					}
					checkSpokeBoundary(-23);
		break;

		case 12: //edge AB
				jpX = parXDash;                                                             
				jpY = componentHeight[compI]*0.5;                                                          
				jpZ = node2Z;                                                                          
				gap = sqrt(pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Edge AB GAP "<<gap<<endl;
					}
					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					boundary = -12;   
					torqueRadius = pXY; //torque due to shear force
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;   
					return isContacted;
				}
				checkSpokeBoundary(-12);
		break;	

		case 26: //edge BF
				jpX = componentWidth[compI]*0.5;                                                             
				jpY = componentHeight[compI]*0.5;                                                          
				jpZ = parZDash;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Edge BF GAP "<<gap<<endl;
					}
					//system("PAUSE");
					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
					jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
					torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force

					boundary = -26;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					torqueRadius = pXY; //torque due to shear force
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;  
					return isContacted;
				}
				checkSpokeBoundary(-26);

		break;
		case 56: //edge EF
				jpX = parXDash;                                                             
				jpY = componentHeight[compI]*0.5;                                                          
				jpZ = node1Z;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				if (gap < 0 )
				{
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Edge EF GAP "<<gap<<endl;
					}
					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);  
					torqueRadius = pXY; //torque due to shear force
					boundary = -56;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;
					return isContacted;
				}
				checkSpokeBoundary(-56);
		break;

		case 34: // edge DC
				jpX = parXDash;                                                             
				jpY = -componentHeight[compI]*0.5;                                                          
				jpZ = node2Z;                                                                          
				gap = sqrt(pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					//system("PAUSE");
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Edge DC GAP "<<gap<<endl;
					}

					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					torqueRadius = pXY; //torque due to shear force
					boundary = -34;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;  
					return isContacted;
				}
				checkSpokeBoundary(-34);
		break;
		case 37: // edge GC
				//jpX = componentWidth[compI]*0.5;                                                             
				//jpY = componentHeight[compI]*0.5;                                                          
				//jpZ = parZDash;                                                                          
				//gap = sqrt(pow(parXDash-jpX,2.0)+pow(parYDash-jpY,2.0)) - rStar; 

				jpX = componentWidth[compI]*0.5;                                                             
				jpY = -componentHeight[compI]*0.5;                                                          
				jpZ = parZDash;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)) - rStar; 
				//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
				//cout<<" GAP "<<gap<<endl;
				if (gap < 0 )
				{
					//system("PAUSE");
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Edge GC GAP "<<gap<<endl;
					}
					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);  
					jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
					jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
					torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force
					boundary = -37;                                                                                                                                                                                                     
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;   
					return isContacted;
				}
				checkSpokeBoundary(-37);
		break;

		case 67: //edge GF
				jpX = componentWidth[compI]*0.5;                                                             
					jpY = parYDash;                                                          
					jpZ = node1Z;                                                                          
					gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
					//cout<<" PAR "<<parXDash<<"   "<<parYDash<<" jpx   jpy  "<<jpX<<"   "<<jpY<<endl;
					//cout<<" GAP "<<gap<<endl;
					if (gap < 0 )
					{
						if (gap < -gapLimit)
						{
							errorFile<<" SPOKE "<<compI<<" Edge GF GAP "<<gap<<endl;
						}
						ijVecX = parXDash - jpX;
						ijVecY = parYDash - jpY;
						ijVecZ = parZDash - jpZ;
						ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

						unitVec[0] = ijVecX/ijXYZ;
						unitVec[1] = ijVecY/ijXYZ;
						unitVec[2] = ijVecZ/ijXYZ;
										//unit vectors for the particle
						ipRVec[0] = -0.5*parDia*(*unitVec);
						ipRVec[1] = -0.5*parDia*(*(unitVec+1));
						ipRVec[2] = -0.5*parDia*(*(unitVec+2));
										//unit vectors for the hole contact point
						jpRVec[0] = jpX;     
						jpRVec[1] = jpY; 
						jpRVec[2] = 0.0; 
																																																														   
						// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
						r1 = ipRVec[0];                                                                                                                                                                                                    
						r2 = ipRVec[1];                                                                                                                                                                                                    
						//r3 = ipRVec[2];                                                                                                                                                                                                  
						ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
						ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
						//ipRVec[2] = does not change                                                                                                                                                                                      
																																																														   
						r1 = jpRVec[0];                                                                                                                                                                                                    
						r2 = jpRVec[1];                                                                                                                                                                                                    
						//r3 = jpRVec[2];                                                                                                                                                                                                  
						jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                        
						jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                        
						//jpRVec[2] = does not change                                                                                                                                                                                      
																																																														   
						r1 = unitVec[0];                                                                                                                                                                                                   
						r2 = unitVec[1];                                                                                                                                                                                                   
						//r3 = unitVec[2];                                                                                                                                                                                                 
						unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
						unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
						jpX = jpX*cos(theta) - jpY*sin(theta)+ ro*cos(theta);
						jpY = jpY*cos(theta) + jpX*sin(theta)+ ro*sin(theta);
						torqueRadius = sqrt(pow(jpX,2)+pow(jpY,2)); //torque due to shear force

						boundary = -67;                                                                                                                                                                                                     
						shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
						findContactForce(1); 
						/*
						if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
						{                                                                                                                                                                                                                  
								checkBoundary(-2);                                                                                                                                                                                             
						}*/                                                                                                                                                                                                                 
						isContacted = 1;  
						return isContacted;
				}
				checkSpokeBoundary(-67);
			break;

			case 78: //edge HG
				jpX = parXDash;                                                             
				jpY = -componentHeight[compI]*0.5;                                                          
				jpZ = node1Z;                                                                          
				gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar; 
				if (gap < 0 )
				{
					//system("PAUSE");
					if (gap < -gapLimit)
					{
						errorFile<<" SPOKE "<<compI<<" Edge HG GAP "<<gap<<endl;
					}

					ijVecX = parXDash - jpX;
					ijVecY = parYDash - jpY;
					ijVecZ = parZDash - jpZ;
					ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

					unitVec[0] = ijVecX/ijXYZ;
					unitVec[1] = ijVecY/ijXYZ;
					unitVec[2] = ijVecZ/ijXYZ;
									//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
									//unit vectors for the hole contact point
					jpRVec[0] = jpX;     
					jpRVec[1] = jpY; 
					jpRVec[2] = 0.0; 
																																																													   
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system                                                                                                                            
					r1 = ipRVec[0];                                                                                                                                                                                                    
					r2 = ipRVec[1];                                                                                                                                                                                                    
					//r3 = ipRVec[2];                                                                                                                                                                                                  
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                         
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                         
					//ipRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = jpRVec[0];                                                                                                                                                                                                    
					r2 = jpRVec[1];                                                                                                                                                                                                    
					//r3 = jpRVec[2];                                                                                                                                                                                                  
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ ro*cos(theta);                                                                                                                                                                         
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ ro*sin(theta);                                                                                                                                                                         
					//jpRVec[2] = does not change                                                                                                                                                                                      
																																																													   
					r1 = unitVec[0];                                                                                                                                                                                                   
					r2 = unitVec[1];                                                                                                                                                                                                   
					//r3 = unitVec[2];                                                                                                                                                                                                 
					unitVec[0] = r1*cos(theta) - r2*sin(theta);                                                                                                                                                                        
					unitVec[1] = r2*cos(theta) + r1*sin(theta);                                                                                                                                                                        
					boundary = -78;
					torqueRadius = pXY;
					shaftAngVel[2] = componentAngVelocity[compI];                                                                                                                                                                      
					findContactForce(1); 
					/*
					if (particleData[parIndex]->SHAFT)                                                                                                                                                                                 
					{                                                                                                                                                                                                                  
							checkBoundary(-2);                                                                                                                                                                                             
					}*/                                                                                                                                                                                                                 
					isContacted = 1;  
					return isContacted;
				}
				checkSpokeBoundary(-78);
				break;
	}
	/*
	if (isContacted == 1)
	{
		nearestSpokeIndex = compI;
	}*/
    return isContacted;
 };

 int Mill::findHubContact(int compI)
 {
	//cout<<" DISK "<<endl;
	calculateTorque = true; //torque calculation is required 
	contact = 0;
    rStar = 0.5*parDia;
	isContacted = 0; //reset
	sfCoff = id->sfricpDisk;
	float radius = 0.0;
	float ijVecX = 0.0;
	float ijVecY = 0.0;
	float ijVecZ = 0.0;
	float ijXYZ  = 0.0;
	float compLength = 0.0;
	
	//float node1X = 0.0; 
	//float node1Y = 0.0;
	//float node1Z = 0.0;
	//float node2X = 0.0;
	//float node2Y = 0.0;
	//float node2Z = 0.0;

	//float theta = 0.0;
	//float ro = 0.0;

	if (fabs(parZ - componentZ[compI]) <= (parDia+componentLength[compI]*0.5)) //particle is close enough to make a contact with disk
	{
		radius = 0.5*componentDia[compI];
		//float theta = componentTheta[compI];
		//float ro = componentR[compI];
		compLength = componentLength[compI];

		// ******** lets take node 1 be the origin of x'y'z' coordinate system
		//float node1X = 0.0; 
		//float node1Y = 0.0;
		float node1Z = 0.0;
		//float node2X = 0.0;
		//float node2Y = 0.0;
		float node2Z = componentLength[compI];
		// ******* new coordinates of the particle in x'y'z' coordinate system
		//parXDash = parX*cos(theta)+parY*sin(theta)- ro;
		//parYDash = parY*cos(theta)-parX*sin(theta);
		parXDash = parX;
		parYDash = parY;
		parZDash = parZ - (componentZ[compI]-compLength*0.5);

		//cout<< "x  "<<xDash<<"   y   "<<yDash<<"   z   "<<zDash<<endl;
		pXY = sqrt(pow(parXDash,2) + pow(parYDash,2)); //Distance to the particle center from z' axis
		if (parZDash > node2Z && parZDash < node2Z+parDia)
		{
			if (pXY < radius)
			{
				contact = 2;
								//errorFile<<"HUB_SIDE_FACE"<<endl;
								//HUB_SIDE_FACE
			}
			else
			{
				contact = 3;
								//errorFile<<"DISK_EDGE_2"<<endl;
								//HUB_EDGE
			}
		}
		else if (parZDash < node1Z && parZDash > node1Z-parDia)
		{
			if (pXY < radius)
			{
				contact = 1;
								//errorFile<<"CYLINDER_SIDE_FACE"<<endl;
								//HUB_SIDE_FACE
			}
			else
			{
				contact = 4;
								//errorFile<<"HUB_EDGE_1"<<endl;
								//CYLINDER_EDGE
			}
		}
		else if (parZDash <= node2Z && parZDash >=  node1Z)
		{
						//cout<<"  axisNode_1[2]   "<<axisNode_1[2]<<endl;
						//cout<<"  axisNode_2[2]   "<<axisNode_2[2]<<endl;
			if (pXY > radius)
			{
				contact = 5;
								//errorFile<<"HUB_OUTERFACE"<<endl;
								//HUB_OUTERFACE
			}
		}

		switch(contact)
		{
			 case 1: //HUB_SIDE_FACE
						  gap = node1Z - parZDash - parDia*0.5;
						  if(gap < 0)
						  {
			  				    if (gap < -gapLimit)
							    {
									errorFile<<" HUB_SIDE_FACE "<<gap<<endl;
									//cout<<" HUB_SIDE_FACE "<<gap<<endl;
								}
							   //cout<<" HUB_SIDE_FACE Overlap  "<<overlap<<endl;
								unitVec[0] = 0.0;
								unitVec[1] = 0.0;
								unitVec[2] = -1.0;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								jpRVec[0] = parXDash;
								jpRVec[1] = parYDash;
								jpRVec[2] = 0.0;
                                boundary = -88;
								torqueRadius = pXY;
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
								if (SHAFT[parIndex])
								{
									checkShaftBoundary(-22);
								}

								isContacted = 1;
								return isContacted;
							}
							checkHubBoundary(-88);

						   break;
						   
						   case 2: //HUB_SIDE_FACE
						   gap = parZDash-node2Z - parDia*0.5;
						   if(gap < 0)
						   {
			  				    if (gap < -gapLimit)
							    {
									errorFile<<" HUB_SIDE_FACE "<<gap<<endl;
									//cout<<" HUB_SIDE_FACE "<<gap<<endl;
								}
								unitVec[0] = 0.0;
								unitVec[1] = 0.0;
								unitVec[2] = 1.0;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								jpRVec[0] = parXDash;
								jpRVec[1] = parYDash;
								jpRVec[2] = 0.0;
								isContacted = 1;
								boundary = -88;
								torqueRadius = pXY;
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
								if (SHAFT[parIndex])
								{
									checkShaftBoundary(-22);
								}
								return isContacted;
								//rotatingSpeed = rotSpeed;
								//cout<<"ipRVec[0]    "<<ipRVec[0]<<endl;
								//cout<<"ipRVec[1]    "<<ipRVec[1]<<endl;
								//cout<<"ipRVec[2]    "<<ipRVec[2]<<endl;     
						   }
						   checkHubBoundary(-88);
						   break;
						   
						   case 3: //HUB_EDGE
						   jpX = parXDash*radius/pXY;
						   jpY = parYDash*radius/pXY;
						   jpZ = node2Z;

						   gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar;
						   if(gap < 0)
						   {
				  				 if (gap < -gapLimit)
							    {
									errorFile<<" HUB_EDGE "<<gap<<endl;
									//cout<<" HUB_EDGE "<<gap<<endl;
								}
								ijVecX = parXDash - jpX;
							   ijVecY = parYDash - jpY;
							   ijVecZ = parZDash - jpZ;
							   ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

								unitVec[0] = ijVecX/ijXYZ;
								unitVec[1] = ijVecY/ijXYZ;
								unitVec[2] = ijVecZ/ijXYZ;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								//unit vectors for the hole contact point
								jpRVec[0] = jpX;     
								jpRVec[1] = jpY; 
								jpRVec[2] = 0.0; 
								
								torqueRadius = radius;
                                boundary = -55;
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);

								isContacted = 1;
								return isContacted;
							}
							checkHubBoundary(-55);
							
						   break;
						   
						   case 4: //HUB_EDGE
						   jpX = parXDash*radius/pXY;
						   jpY = parYDash*radius/pXY;
						   jpZ = node1Z;

						   gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar;
						   if(gap < 0)
						   {
				  				if (gap < -gapLimit)
							    {
									errorFile<<" HUB_EDGE "<<gap<<endl;
									//cout<<" HUB_EDGE "<<gap<<endl;
								}
								ijVecX = parXDash - jpX;
							   ijVecY = parYDash - jpY;
							   ijVecZ = parZDash - jpZ;
							   ijXYZ = sqrt(pow(ijVecX,2)+pow(ijVecY,2)+pow(ijVecZ,2));

								unitVec[0] = ijVecX/ijXYZ;
								unitVec[1] = ijVecY/ijXYZ;
								unitVec[2] = ijVecZ/ijXYZ;
								//cout<<" unitVec  "<<unitVec[2]<<endl;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								jpRVec[0] = jpX;     
								jpRVec[1] = jpY; 
								jpRVec[2] = 0.0; 
								
								torqueRadius = radius;
                                boundary = -55;
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
							
								isContacted = 1;
								return isContacted;
						   }
						   checkHubBoundary(-55);
						   break;
						   
						   case 5: //HUB_OUTERFACE
						   gap = pXY - radius - rStar;
						   if(gap < 0)
						   {
				  				if (gap < -gapLimit)
							    {
									errorFile<<" HUB_OUTERFACE "<<gap<<endl;
									//cout<<" HUB_OUTERFACE "<<gap<<endl;
								}
								unitVec[0] = parXDash/pXY;
								unitVec[1] = parYDash/pXY;
								unitVec[2] = 0.0;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								//cout<< "ipRVec[1]  "<<ipRVec[1]<<endl;
								jpRVec[0] = radius*(*unitVec);     
								jpRVec[1] = radius*(*(unitVec+1)); 
								jpRVec[2] = radius*(*(unitVec+2)); 
                                //errorFile<<" U VEC "<<unitVec[0]<<"   "<<unitVec[1]<<endl;
								torqueRadius = radius;
								boundary = -44;
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
								isContacted = 1;
								return isContacted;
						   }
						   checkHubBoundary(-44);
						   //cout<<" gap  "<<gap<<endl;
						  
					 break;
		}//end of switch
	}//particle close enough
	return isContacted;
 };

int Mill::findShaftContact(int compI)
{
	isContacted = 0; //reset
	rStar = 0.5*parDia;
	pXY = sqrt(pow(parX,2)+pow(parY,2));
	gap = pXY - 0.5*componentDia[compI] - rStar;
	sfCoff = id->sfricpDisk;
	if(gap < 0)
	{
		if (gap < -gapLimit)
		{
			errorFile<<"gap "<<gap<<" SHAFT"<<endl;
			cout<<"gap "<<gap<<" SHAFT"<<endl;
		}
		unitVec[0] = parX/pXY;//always towards the centre from the contact point//
		unitVec[1] = parY/pXY;
		unitVec[2] = 0.0;

		//File<<" UNIT "<<unitVec[0]<<"  "<<unitVec[1]<<"  "<<unitVec[2];
		ipRVec[0] = -0.5*parDia*(*unitVec);
		ipRVec[1] = -0.5*parDia*(*(unitVec+1));
		ipRVec[2] = -0.5*parDia*(*(unitVec+2));

		jpRVec[0] = 0.5*componentDia[compI]*(*unitVec);     
		jpRVec[1] = 0.5*componentDia[compI]*(*(unitVec+1)); 
		jpRVec[2] = 0.5*componentDia[compI]*(*(unitVec+2)); 
		boundary = -2;
		torqueRadius = 0.5*componentDia[compI]; //distance to the contact point
		//compType = 2; //for torque calculation  
	
		shaftAngVel[2] = componentAngVelocity[compI];
	    findContactForce(1);
		if (DRUM_SIDE_FACE[parIndex])
		{
			checkDrumBoundary(-33);
		}

		if (HUB_SIDE_FACE[parIndex])
		{
			checkHubBoundary(-88);
		}
		isContacted = 1;
		return isContacted;
	}
	checkShaftBoundary(-22); //check for shaft itself
	return isContacted;
};

int Mill::findDrumContact(int compI)
{
	isContacted = 0; //reset
	rStar = 0.5*parDia;
	pXY = sqrt(pow(parX,2)+pow(parY,2));
	gap = 0.5*componentDia[compI] - pXY - rStar;
	sfCoff = id->sfricpWall; 
	//errorFile<<" DRUM "<<endl;
	//File<<" GAP "<<gap;
	if(gap < 0)
	{
        if (gap < -gapLimit)
		{
			errorFile<<"gap "<<gap<<" DRUM"<<endl;
			cout<<"gap "<<gap<<" DRUM"<<endl;
		}
		unitVec[0] = -parX/pXY;//always towards the centre from the contact point//
		unitVec[1] = -parY/pXY;
		unitVec[2] = 0.0;

		ipRVec[0] = -0.5*parDia*(*unitVec);
		ipRVec[1] = -0.5*parDia*(*(unitVec+1));
		ipRVec[2] = -0.5*parDia*(*(unitVec+2));

		jpRVec[0] = -0.5*componentDia[compI]*(*unitVec);
		jpRVec[1] = -0.5*componentDia[compI]*(*(unitVec+1));
		jpRVec[2] = -0.5*componentDia[compI]*(*(unitVec+2));
		boundary = -11;
		shaftAngVel[2] = componentAngVelocity[compI];
		findContactForce(1);
        //errorFile<<" DRUM "<<shaftAngVel[2]<<endl;
		if (DRUM_SIDE_FACE[parIndex])
		{
			checkDrumBoundary(-33); //check for end wall
		}

		isContacted = 1;
		return isContacted;
	}
	checkDrumBoundary(-11); //check for drum
	return isContacted;
};

int Mill::findFirstEndWallContact(int compI)
{
	isContacted = 0;
	sfCoff = id->sfricpWall; 
	if (parZ - (componentZ[compI]+componentLength[compI]*0.5) < parDia/2.0)
	{	
		//case 1://CONTACT_WITH_DRUM_SIDEFACE
		unitVec[0] = 0.0;
		unitVec[1] = 0.0;
		unitVec[2] = 1.0;
		gap = parZ - parDia/2.0;
		if(gap < 0)
		{	
           if (gap < -gapLimit)
			{
				errorFile<<"gap "<<gap<<" DRUM SIDEFACE"<<endl;
				cout<<"gap "<<gap<<" DRUM SIDEFACE"<<endl;
			}
			//contactString = "CONTACT_WITH_DRUM_SIDEFACE";
			ipRVec[0] = -0.5*parDia*(*unitVec);
			ipRVec[1] = -0.5*parDia*(*(unitVec+1));
			ipRVec[2] = -0.5*parDia*(*(unitVec+2));
			jpRVec[0] = parX;
			jpRVec[1] = parY;
			jpRVec[2] = 0.0;
			boundary = -33;
			shaftAngVel[2] = componentAngVelocity[compI];
			findContactForce(1);
			if (DRUM_OUTER_FACE[parIndex])
			{
				checkDrumBoundary(-11);
			}
			if (SHAFT[parIndex])
			{
				checkShaftBoundary(-22);
			}
			isContacted = 1;
			return isContacted;
		}
	}
	checkDrumBoundary(-33);
	return isContacted;
};

int Mill::findSecondEndWallContact(int compI)
{
	isContacted = 0;
	sfCoff = id->sfricpWall; 
	if (parZ > (componentZ[compI]-componentLength[compI]*0.5 - parDia/2.0) )
	{	
		//case 1://CONTACT_WITH_DRUM_SIDEFACE
		unitVec[0] = 0.0;
		unitVec[1] = 0.0;
	
		unitVec[2] = -1.0;
		gap = (componentZ[compI]-componentLength[compI]*0.5) - parZ - parDia/2.0 ;
		if(gap < 0)
		{	
			if (gap < -gapLimit)
			{
				errorFile<<"gap "<<gap<<" DRUM SIDEFACE"<<endl;
				cout<<"gap "<<gap<<" DRUM SIDEFACE"<<endl;
			}
			//contactString = "CONTACT_WITH_DRUM_SIDEFACE";
				
			ipRVec[0] = -0.5*parDia*(*unitVec);
			ipRVec[1] = -0.5*parDia*(*(unitVec+1));
			ipRVec[2] = -0.5*parDia*(*(unitVec+2));
			jpRVec[0] = parX;
			jpRVec[1] = parY;
			jpRVec[2] = 0.0;
			boundary = -33;
			shaftAngVel[2] = componentAngVelocity[compI];
			findContactForce(1);
			if (DRUM_OUTER_FACE[parIndex])
			{
				checkDrumBoundary(-11);
			}
			if (SHAFT[parIndex])
			{
				checkShaftBoundary(-22);
			}
			isContacted = 1;
			return isContacted;
		}
	}
	checkDrumBoundary(-33);
	return isContacted;
};
//-------------------------------------------------------------------

bool Mill::checkData(int showKE, int parI, int arrI, float x,float y,float z,float cF, 
				float sF,float sE, float gap, float kE,float rVel, float tW)
{
	bool found = false; // To check contact is existing or not
	for (int i=0; i<noOfContacts[parI] ; i++)
	{
		if (contactList[parI][i] == arrI) // Contact  exists
		{	
			found = true;
		}
	}
		
	if (!found) // A new contact
	{ 
		switch(arrI)
		{
			
			case -11:
				DRUM_OUTER_FACE[parI] = true;
				break;
			case -22:
				SHAFT[parI] = true;
				break;
			case -33:
				DRUM_SIDE_FACE[parI] = true;
				break;
			case -44:
				HUB_OUTER_FACE[parI] = true;
				break;
			case -88:
				HUB_SIDE_FACE[parI] = true;
				break;
		}
		 totalWear[parI] = tW;
		contactList[parI][noOfContacts[parI]] = arrI;
		noOfContacts[parI] = noOfContacts[parI] + 1;
		totalCollisions[parI] = totalCollisions[parI] + 1;//keeps a track of total number of collisions, particle having up until any time
		if(showKE == 1)
		{
			forceFile.write(reinterpret_cast <const char*>(&kE),sizeof(float));
		}
	}
	return found;
}


bool Mill::deleteContact(int parI, int arrIndex)
{
	bool deleted = false;
	for (int i=0; i<noOfContacts[parI]; i++ )
	{
		if (contactList[parI][i] == arrIndex)
		{
			contactList[parI][i] = contactList[parI][noOfContacts[parI] -1];
			contactList[parI][noOfContacts[parI] -1] = -10;  
			noOfContacts[parI] = noOfContacts[parI] - 1;
			deleted = true;
			break;
			
		}
	}
	
	switch(arrIndex)
	{
			case -11:
				DRUM_OUTER_FACE[parI] = false;
				break;
			case -22:
				SHAFT[parI] = false;
				break;
			case -33:
				DRUM_SIDE_FACE[parI] = false;
				break;
			case -44:
				HUB_OUTER_FACE[parI] = false;
				break;
			case -88:
				HUB_SIDE_FACE[parI] = false;
				break;
	}
	return deleted;
};


void Mill::checkHubBoundary(int prevB)
{
    rStar = 0.5*particleDiameter[parIndex];
	pXY = sqrt(pow(parX,2)+pow(parY,2));
	gap = 0.5*id->chamberInnerDia - pXY - rStar;
    switch(prevB)
	{
		case -44: // hub outer face
		gap = pXY - 0.5*hubDia - rStar;
        if (gap > collisionGapCoeff)
		{
			deleteContact(parIndex, prevB);
		}
		break;

		case -55: // hub outer line
		jpX = (parX*hubDia/2.0)/pXY;//X coordinate on disk outline//
		
		jpY = (parY*hubDia/2.0)/pXY;//Y coordinate on disk outline//
		if(parZ - nearestHubZ > 0)
			jpZ = hubThick/2.0;//Z coordinate on disk outline//
		else
			jpZ = -hubThick/2.0;//Z coordinate on disk outline//
		ijVec[0] = parX-jpX;
		ijVec[1] = parY-jpY;
		ijVec[2] = parZ - nearestHubZ - jpZ;
		unitVector(ijVec);
		unitVec[0] = r1;
		unitVec[1] = r2;
		unitVec[2] = r3;
		
		gap = sqrt(pow(parX-jpX,2)+pow(parY-jpY,2)+pow(parZ - nearestHubZ - jpZ,2)) - rStar;
		if (gap > collisionGapCoeff)
		{
			deleteContact(parIndex, prevB);
		}
		break;

		case -88: // hub side face
		unitVec[0] = 0.0;
		unitVec[1] = 0.0;
		if(parZ > nearestHubZ)
		{
			unitVec[2] =1.0;
		}	
		else
		{
			unitVec[2] = -1.0;
		}
			
		gap = fabs(parZ - nearestHubZ) - (hubThick/2.0+parDia/2.0);
		if (gap > collisionGapCoeff)
		{
			deleteContact(parIndex, prevB);
		}
		break;
	}
	
};


void Mill::checkDrumBoundary(int prevB)
{
	rStar = 0.5*particleDiameter[parIndex];
	pXY = sqrt(pow(parX,2)+pow(parY,2));
	gap = 0.5*id->chamberInnerDia - pXY - rStar;

	switch(prevB)
	{
		case -11: // drum
		if (gap > collisionGapCoeff)
		{
			deleteContact(parIndex, prevB);
		}
		break;
		case -33: // drum side face
		unitVec[0] = 0.0;
		unitVec[1] = 0.0;
		
		if(parZ - parDia/2.0 > 0)
		{
			unitVec[2] = -1.0;
			gap = id->chamberLength - parZ - parDia/2.0 ;
		}
		else
		{
			unitVec[2] = 1.0;
			gap = parZ - parDia/2.0;
		}
		if (gap > collisionGapCoeff)
		{
			deleteContact(parIndex, prevB);
		}
		break;
	}

};

void Mill::checkShaftBoundary(int prevB)
{
	rStar = 0.5*particleDiameter[parIndex];
	pXY = sqrt(pow(parX,2)+pow(parY,2));
	gap = pXY - 0.5*shaftDia - rStar;
	if (gap > collisionGapCoeff)
	{
		deleteContact(parIndex, prevB);
	}
};

void Mill::checkSpokeBoundary(int prevB)
{
	ipDia = particleDiameter[parIndex];
	rStar = 0.5*particleDiameter[parIndex];
    pXY = sqrt(pow(parX,2)+pow(parY,2));
	float node1Z = 0.0;
	float node2Z = componentLength[nearestSpokeIndex];
	float theta = componentTheta[nearestSpokeIndex];
	float ro = componentR[nearestSpokeIndex];
	float compLength = componentLength[nearestSpokeIndex];
	float compHeight = componentHeight[nearestSpokeIndex];
	float compWidth = componentWidth[nearestSpokeIndex];
	// ******* new coordinates of the particle in x'y'z' coordinate system
	parXDash = parX*cos(theta)+parY*sin(theta)- ro;
	parYDash = parY*cos(theta)-parX*sin(theta);
	parZDash = parZ - (componentZ[nearestSpokeIndex]-compLength*0.5);
	
	switch(prevB)
	{
		case -1256: //top plane
			if (parYDash - parDia*0.5 - compHeight*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);
			}
		break;
	
		case -3478: //bottom plane
			if (fabs(parYDash) - parDia*0.5 - compHeight*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);
			}
		break;

		case -2367: //plane BCGF
			if (parXDash - parDia*0.5 - compWidth*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -1234: //plane ABCD
			if (parZDash - parDia*0.5 - compLength > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;


		case -5678: //plane EFGH
			if (fabs(parZDash) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;
		
		case -12: //edge AB
			if (sqrt(pow((node2Z-parZDash),2)+pow((compHeight*0.5-parYDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -56: //edge EF
			if (sqrt(pow((node1Z-fabs(parZDash)),2)+pow((compHeight*0.5-parYDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -34: //edge CD
			if (sqrt(pow((node2Z-parZDash),2)+pow((compHeight*0.5-fabs(parYDash)),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -78: //edge GH
			if (sqrt(pow((node1Z-fabs(parZDash)),2)+pow((compHeight*0.5-fabs(parYDash)),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;


		case -23: //edge BC
			if (sqrt(pow((node2Z-parZDash),2)+pow((compWidth*0.5-parXDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -67: //edge FG
			if (sqrt(pow((node1Z-fabs(parZDash)),2)+pow((compWidth*0.5-parXDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -26: //edge BF
			if (sqrt(pow((compHeight*0.5-parYDash),2)+pow((compWidth*0.5-parXDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -37: //edge CG
			if (sqrt(pow((compHeight*0.5-fabs(parYDash)),2)+pow((compWidth*0.5-parXDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -2: //corner B
			if (sqrt(pow((compHeight*0.5-parYDash),2)+pow((compWidth*0.5-parXDash),2)+pow((node2Z-parZDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -3: //corner C
			if (sqrt(pow((compHeight*0.5-fabs(parYDash)),2)+pow((compWidth*0.5-parXDash),2)+pow((node2Z-parZDash),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -6: //corner F
			if (sqrt(pow((compHeight*0.5-parYDash),2)+pow((compWidth*0.5-parXDash),2)+pow((node1Z-fabs(parZDash)),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;

		case -7: //corner G
			if (sqrt(pow((compHeight*0.5-fabs(parYDash)),2)+pow((compWidth*0.5-parXDash),2)+pow((node1Z-fabs(parZDash)),2)) - parDia*0.5 > collisionGapCoeff)
			{
				deleteContact(parIndex, prevB);	
			}
		break;
	
	}
};

void Mill::findContactForce(int contact)
{
	bool contacted = false;
	if (gap < 0 )
	{
		contacted = true; //a contact exists
		float velArray[3];
		velArray[0] = particleAngVelX[parIndex];
		velArray[1] = particleAngVelY[parIndex];
		velArray[2] = particleAngVelZ[parIndex];
		crossProduct(velArray,ipRVec);
		rotVel[0] = r1;
		rotVel[1] = r2;
		rotVel[2] = r3;

		velArray[0] = particleVelX[parIndex];
		velArray[1] = particleVelY[parIndex];
		velArray[2] = particleVelZ[parIndex];
		nrmVel = dotProduct(velArray,unitVec);
		
		vecAdd(velArray,rotVel);
		ipCntPntVel[0] = r1;     
		ipCntPntVel[1] = r2; 
		ipCntPntVel[2] = r3; 

		crossProduct(shaftAngVel,jpRVec);//shaft and disks rotate at the same speed//
		jpCntPntVel[0] = r1;  
		jpCntPntVel[1] = r2;  
		jpCntPntVel[2] = r3;  

		vecSubstract(ipCntPntVel,jpCntPntVel);
		cntPntVel[0] = r1;  
		cntPntVel[1] = r2;  
		cntPntVel[2] = r3;
		instRelativeVel = dotProduct(cntPntVel,unitVec);
		normDampC = id->pdNormDampC;
		calculateForce();
		//---Torque calculation is only for the components 
		//--- which are mounted on the shaft 
		if(nrmForce < 0)
        {
            nrmForce = 0.0;
        }
        float forceX = -unitVec[0]*nrmForce;
		float forceY = -unitVec[1]*nrmForce;
		float forceZ = -unitVec[2]*nrmForce;

        float resultantF = sqrt(pow(forceX,2)+pow(forceY,2));
       
		switch(boundary)
		{
			case -22: //shaft
			torque = torque + resultantF*torqueRadius*sfCoff;
			break;
			
			case -44: //disk outer face
			torque = torque + resultantF*torqueRadius*sfCoff;
			break;
			
			case -55: //disk edge
			torque = torque + resultantF*torqueRadius*sfCoff;
			break;
			
			case -88: //disk side face
			torque = torque + fabs(forceZ)*torqueRadius*sfCoff;
			break;
			//------------------ Pin Contacts --------------------
			case -1234: //plane ABCD (side face)
			torque = torque + fabs(forceZ)*torqueRadius*sfCoff;
			break;

			case -5678: //plane EFGH (side face)
			torque = torque + fabs(forceZ)*torqueRadius*sfCoff;
			break;

			default: //for all other contacts 
			torque = torque - forceX*collPointY + forceY*collPointX;
		}

        /*
        if (compType == 1 || compType == 2) //disks, holes and shaft lies in between these type
		{
			crossProduct(jpRVec,rotorForce);
			momentum[0] = r1;  
			momentum[1] = r2;  
			momentum[2] = r3;  
			
			temp = sfCoff*nrmCntForce*torqueRadius;
			momentum[2] = *(momentum+2) + temp;
			//errorFile<<" mom  "<<jpRVec[0]<<"   "<<jpRVec[1]<<"   "<<jpRVec[2]<<endl;		
			vecAdd(momentum,torqueOnMill);
			torqueOnMill[0] = r1;  
			torqueOnMill[1] = r2;  
			torqueOnMill[2] = r3;  

		}*/
    }
	//---- End of torque calculation
	if (!contacted)
	{
		//contactString = "CONTACT_WITH_NONE";
		disp[0] = 0.0;
		disp[1] = 0.0;
		disp[2] = 0.0;
		particleHisDispX[parIndex] = 0.0;
		particleHisDispY[parIndex] = 0.0;
		particleHisDispZ[parIndex] = 0.0;
		//checkAllBoundary();
	}
 };

