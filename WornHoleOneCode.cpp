#include "Mill.h"

void Mill::generateDisc()
{
	int hollowCompIndex = 0; //array index for number of hollow components 

	for (int i=0; i<id->noOfDisks ;i++ )
	{
		hollowCompIndex = 0;
		for (int j=0; j<id->noOfHoles ;j++ )
		{
			hollowCompR[compIndex][hollowCompIndex] = holeCenterRadius*id->lengthFactor;
			hollowCompTheta[compIndex][hollowCompIndex] = 2.0*M_PI*j/id->noOfHoles - M_PI/id->noOfHoles;
			//hollowCompZ[compIndex][hollowCompIndex] = (id->diskApart+id->diskThick)*(0.5+i);//-- WIth boundary conditions ------
			hollowCompZ[compIndex][hollowCompIndex] = 0.5*(id->diskApart+id->diskThick) + i*(id->diskApart+id->diskThick);
			hollowCompLength[compIndex][hollowCompIndex] = id->diskThick;
			hollowCompDia[compIndex][hollowCompIndex] = id->holeDia;
			//subComponentType[compIndex][hollowCompIndex] = 1; //type is hollow cylinder
			//hollowCompAngVelocity[compIndex][hollowCompIndex] = 0.0; //mill rotation speed
			hollowCompAngVelocity[compIndex][hollowCompIndex] = -2.0*M_PI*id->rotSpeed/60.0; //mill rotation speed
			hollowCompIndex++;
		}
        
	
	noOfSubComponents[compIndex] = hollowCompIndex;
        componentType[compIndex] = 1; //type is disk
        componentR[compIndex] = 0.0;
        componentDia[compIndex] = id->diskDia;
        componentTheta[compIndex] = 0.0;
        //componentZ[compIndex] = (id->diskApart+id->diskThick)*(0.5+i);//-- WIth boundary conditions ------
        componentZ[compIndex] = 0.5*(id->diskApart+id->diskThick) + i*(id->diskApart+id->diskThick);
        componentLength[compIndex] = id->diskThick;
        componentAngVelocity[compIndex] = -2.0*M_PI*id->rotSpeed/60.0; //mill rotation speed
        compIndex++;
	}
	//------------------------------------
};

int Mill::findCompContact(int compI)
{
	int type = componentType[compI];
	compType = 0;
	switch(type)
	{
		case 1: //disk
		findDiskContact(compI);
		break;

		case 2: //shaft
		findShaftContact(compI);
		break;

		case 3: //drum
		findDrumContact(compI);
		break;

		case 4: //first end wall (z=0)
		findFirstEndWallContact(compI);
		break;

		case 5: //second end wall 
		findSecondEndWallContact(compI);
		break;
	
	}
    return 1;
};

int Mill::findDiskContact(int compI)
{
    //cout<<" DISK "<<endl;
	contact = 0;
    rStar = 0.5*parDia;
	isInside = false;  //particle is inside the hole
	isContacted = 0; //reset
	sfCoff = id->sfricpDisk;
	double prevTempXY = 10000.0*id->lengthFactor; //previous XY distance to the vacant region 
	double tempXY = 0.0;
	int nearestHoleIndex = 0;
	double radius = 0.0;
	double ijVecX = 0.0;
	double ijVecY = 0.0;
	double ijVecZ = 0.0;
	double ijXYZ  = 0.0;
	double compLength = 0.0;
	
	double node1X = 0.0; 
	double node1Y = 0.0;
	double node1Z = 0.0;
	double node2X = 0.0;
	double node2Y = 0.0;
	double node2Z = 0.0;

	double theta = 0.0;
	double ro = 0.0;

	if (fabs(parZ - componentZ[compI]) <= (parDia+componentLength[compI]*0.5)) //particle is close enough to make a contact with disk
	{
		for (int j=0; j<noOfSubComponents[compI] ;j++ ) //for each hollow component find the nearest hollow region
		{
			//parXDash, parYDash and parZDash are the new coordinates of the particle 
			//considering the hollow region center as the origin of the new coordinate system
			//parXDash = x*cos(theta)+y*sin(theta)-ro;
			//parYDash = x*sin(theta)+y*sin(theta);
			theta = hollowCompTheta[compI][j];
            ro = hollowCompR[compI][j];
			double holeX = ro*cos(theta);
            double holeY = ro*sin(theta);
            
			tempXY = sqrt(pow((parX - holeX),2)+pow((parY - holeY),2)); //distance between particle center and hole center
			if (tempXY < prevTempXY)
			{
			    prevTempXY = tempXY;
				nearestHoleIndex = j;
				nearestHoleX = holeX;
				nearestHoleY = holeY;
				nearestDiskZ = componentZ[compI]; //nearest disk
			}
		}
        radius = 0.5*hollowCompDia[compI][nearestHoleIndex];
		theta = hollowCompTheta[compI][nearestHoleIndex];
		ro = hollowCompR[compI][nearestHoleIndex];
		compLength = hollowCompLength[compI][nearestHoleIndex];
		//--------- lets take node 1 be the origin of x'y'z' coordinate system
		node1X = 0.0;                               
		node1Y = 0.0;                               
		node1Z = 0.0;                               
		node2X = 0.0;                               
		node2Y = 0.0;                               
		node2Z = compLength; 
		parXDash = parX*cos(theta)+parY*sin(theta)- ro;
		parYDash = parY*cos(theta)-parX*sin(theta);
		parZDash = parZ - (hollowCompZ[compI][nearestHoleIndex]-compLength*0.5);
		pXY = sqrt(pow(parXDash,2)+pow(parYDash,2));
		if(parZDash >= 0 && parZDash <= node2Z) // partilce is inside the hole
		{     
            if(pXY < radius)
            {
                 contact = 3;//CONTACT_WITH_HOLE_FACE
                isInside = true;//particle is inside
            }
        }
		else if(parZDash - node2Z > 0 && parZDash - node2Z < parDia*0.5)
		{
				 if(pXY < radius)
				{
					   contact = 2;//CONTACT_WITH_HOLE_EDGE 2
					   //errorFile<<" HOLE EDGE 2"<<endl;
					   isInside = true;
				}
		}
		else if(node1Z - parZDash > 0 && node1Z - parZDash < parDia*0.5)
		{
				if(pXY < radius)
				{
						contact = 1;//CONTACT_WITH_HOLE_EDGE 1
						//errorFile<<" HOLE EDGE 1"<<endl;
						isInside = true;
				}
		}

		switch(contact)
		{
			   case 1: //CONTACT_WITH_HOLE_EDGE 1
			   xx = parXDash - node1X;
			   yy = parYDash - node1Y;
			   //zz = parZDash - node1Z;
			   zz = 0.0;
			   pXY = sqrt(pow(xx, 2)+pow(yy, 2));
			   if (pXY != 0)
			   {
					unitVec[0] = xx/pXY;
					unitVec[1] = yy/pXY;
					unitVec[2] = zz/pXY;
			   }
			   else
   			   {
					unitVec[0] = 0.0;
					unitVec[1] = -1.0;
					unitVec[2] = 0.0;
			   }
			   
			   jpX = node1X + radius*(*unitVec);
			   jpY = node1Y + radius*(*(unitVec+1));
			   jpZ = node1Z +  radius*(*(unitVec+2));;
			   gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar;
              
			   if (gap < 0)
			   {    
				    if (gap < -gapLimit)
				    {
						errorFile<<"CONTACT_WITH_HOLE_EDGE 1     "<<gap/id->lengthFactor<<endl;
						cout<<"CONTACT_WITH_HOLE_EDGE 1     "<<gap/id->lengthFactor<<endl;
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
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ro*cos(theta);
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ro*sin(theta);
					//jpRVec[2] = does not change
				
					r1 = unitVec[0];
					r2 = unitVec[1];
					//r3 = unitVec[2];
					unitVec[0] = r1*cos(theta) - r2*sin(theta);
					unitVec[1] = r2*cos(theta) + r1*sin(theta);
					//unitVec[2] = does not change
					boundary = -7;
					torqueRadius = pXY;
					compType = 1; //for torque calculation hole is a part of disk 
					shaftAngVel[2] = hollowCompAngVelocity[compI][nearestHoleIndex];
					findContactForce(1);

					isContacted = 1;
					isInside = true;
					return isContacted;
			   }
			   checkBoundary(-7); //check for disk hole line
			   break;
			   
			   case 2: //CONTACT_WITH_HOLE_EDGE 2
			   xx = parXDash - node2X;
			   yy = parYDash - node2Y;
			   zz = 0.0;
			   pXY = sqrt(pow(xx, 2)+pow(yy, 2));
			   if (pXY != 0)
			   {
					unitVec[0] = xx/pXY;
					unitVec[1] = yy/pXY;
					unitVec[2] = zz/pXY;
			   }
			   else
   			   {
					unitVec[0] = 0.0;
					unitVec[1] = -1.0;
					unitVec[2] = 0.0;
			   }

			   jpX = node2X + radius*(*unitVec);
			   jpY = node2Y + radius*(*(unitVec+1));
			   jpZ = node2Z;
			   gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar;
			   if (gap < 0)
			   {
   				    if (gap < -gapLimit)
				    {
				    	errorFile<<"CONTACT_WITH_HOLE_EDGE 2      "<<gap/id->lengthFactor<<endl;
						cout<<"CONTACT_WITH_HOLE_EDGE 2      "<<gap/id->lengthFactor<<endl;
					}
			        //errorFile<<" HOLE EDGE 2"<<"  X  "<<parXDash<<" Y  "<<parYDash<<" Z "<<parZDash<<endl;
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
					jpRVec[0] = r1*cos(theta) - r2*sin(theta)+ro*cos(theta);
					jpRVec[1] = r2*cos(theta) + r1*sin(theta)+ro*sin(theta);
					//jpRVec[2] = does not change
				
					r1 = unitVec[0];
					r2 = unitVec[1];
					//r3 = unitVec[2];
					unitVec[0] = r1*cos(theta) - r2*sin(theta);
					unitVec[1] = r2*cos(theta) + r1*sin(theta);
					//unitVec[2] = does not change
					boundary = -7;
					torqueRadius = pXY;
					compType = 1; //for torque calculation hole is a part of disk 
					shaftAngVel[2] = hollowCompAngVelocity[compI][nearestHoleIndex];
					findContactForce(1);

					isContacted = 1;
					isInside = true;
					return isContacted;
			   }
			   checkBoundary(-7); //check for disk hole line
			   break;
			   
			   case 3: //CONTACT_WITH_HOLE_FACE
			   pXY = sqrt(pow(parXDash,2)+pow(parYDash,2));
               gap = radius - pXY - rStar;
			   
			   if(gap < 0)
			   {
				   //if (pXY == 0.0)
				   //{	
					//   cout<<" PXY "<<pXY<<endl;
				   //}
   				    if (gap < -gapLimit)
				    {
					   errorFile<<"CONTACT_WITH_HOLE_FACE      "<<gap/id->lengthFactor<<endl;
					    cout<<"CONTACT_WITH_HOLE_FACE      "<<gap/id->lengthFactor<<endl;
					}			   
					unitVec[0] = - parXDash/pXY;
					unitVec[1] = - parYDash/pXY;
					unitVec[2] = 0.0;
					
					//unit vectors for the particle
					ipRVec[0] = -0.5*parDia*(*unitVec);
					ipRVec[1] = -0.5*parDia*(*(unitVec+1));
					ipRVec[2] = -0.5*parDia*(*(unitVec+2));
					//unit vectors for the hole contact point
					double contactX = -radius*(*unitVec);
					double contactY = -radius*(*(unitVec+1));
					double contactZ = -radius*(*(unitVec+2));
					jpRVec[0] = contactX; 
					jpRVec[1] = contactY; 
					jpRVec[2] = contactZ; 
					
					// Now convert ipRVec, jpRvec, and unitVector to the original coordinate (x,y,z) system
					r1 = ipRVec[0];
					r2 = ipRVec[1];
					//r3 = ipRVec[2];
					//ipRVec[0] = r1*cos(theta) + r2*sin(theta) + ro*cos(theta);
					//ipRVec[1] = r2*cos(theta) - r1*sin(theta) - ro*sin(theta);
					ipRVec[0] = r1*cos(theta) - r2*sin(theta);
					ipRVec[1] = r2*cos(theta) + r1*sin(theta);

					//ipRVec[2] = does not change
						
					r1 = jpRVec[0];
					r2 = jpRVec[1];
					//r3 = jpRVec[2];
					jpRVec[0] = r1*cos(theta) - r2*sin(theta) + ro*cos(theta) ;
					jpRVec[1] = r2*cos(theta) + r1*sin(theta) + ro*sin(theta);
					//jpRVec[2] = does not change
				
					r1 = unitVec[0];
					r2 = unitVec[1];
					//r3 = unitVec[2];
					unitVec[0] = r1*cos(theta) - r2*sin(theta);
					unitVec[1] = r2*cos(theta) + r1*sin(theta);
					
					boundary = -6;
					torqueRadius = pXY;
					compType = 1; //for torque calculation hole is a part of disk 

					shaftAngVel[2] = hollowCompAngVelocity[compI][nearestHoleIndex];
					findContactForce(1);
					//unitVec[2] = does not change
					isContacted = 1;
					isInside = true;
                    //errorFile<<" HOLE FACE "<<endl;
					return isContacted;
			  }
			   checkBoundary(-6); //check for disk hole face
               break;

		}//end of switch (end of hole contact calculation)
		if (isInside == false) //does not contact with hole
		{
				  	radius = 0.5*componentDia[compI];
					theta = componentTheta[compI];
					ro = componentR[compI];
					compLength = componentLength[compI];

					//-------- lets take node 1 be the origin of x'y'z' coordinate system
					node1X = 0.0; 
					node1Y = 0.0;
					node1Z = 0.0;
					node2X = 0.0;
					node2Y = 0.0;
					node2Z = componentLength[compI];
					// ******* new coordinates of the particle in x'y'z' coordinate system
					//parXDash = parX*cos(theta)+parY*sin(theta)- ro;
					//parYDash = parY*cos(theta)-parX*sin(theta);
					//parZDash = parZ - (componentZ[compI]-compLength*0.5);
					parXDash = parX;
					parYDash = parY;
					parZDash = parZ - (componentZ[compI]-compLength*0.5);

					pXY = sqrt(pow(parXDash,2) + pow(parYDash,2)); //Distance to the particle center from z' axis
					contact = 0;
					if (parZDash > node2Z && parZDash < node2Z+parDia)
					{
						if (pXY < radius)
						{
							contact = 2;
								//CYLINDER_SIDE_FACE
						}
						else
						{
							contact = 3;
									//CYLINDER_EDGE
						}
					}
					else if (parZDash < node1Z && parZDash > node1Z-parDia)
					{
						if (pXY < radius)
						{
							contact = 1;
								//DISK_SIDE_FACE
						}
						else
						{
							contact = 4;
								//CYLINDER_EDGE
						}
					}
					else if (parZDash <= node2Z && parZDash >=  node1Z)
					{
						if (pXY > radius)
						{
							contact = 5;
								//DISK_OUTERFACE
						}
					}

					switch(contact)
					{
						  case 1: //DISK_SIDE_FACE
						  gap = node1Z - parZDash - parDia*0.5;
						  if(gap < 0)
						  {
   								if (gap < -gapLimit)
							    {
									errorFile<<" DISK_SIDE_FACE "<<gap/id->lengthFactor<<endl;
									cout<<" DISK_SIDE_FACE "<<gap/id->lengthFactor<<endl;
								}
							    //errorFile<<" DISK_SIDE_FACE 1 Overlap  "<<endl;
								unitVec[0] = 0.0;
								unitVec[1] = 0.0;
								unitVec[2] = -1.0;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								//unit vectors for the hole contact point
								//contactVector[0] = -radius*(*unitVec);
								//contactVector[1] = -radius*(*(unitVec+1));
								//contactVector[2] = -radius*(*(unitVec+2));
								jpRVec[0] = parXDash;
								jpRVec[1] = parYDash;
								jpRVec[2] = 0.0;
                                boundary = -8;
								torqueRadius = pXY; //distance to the contact point
								compType = 1; //for torque calculation  
			
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
								if (SHAFT[parIndex])
								{
									checkBoundary(-2);
								}

								isContacted = 1;
								return isContacted;
							}
						   checkBoundary(-8); //check for disk side face
						   break;
						   
						   case 2: //DISK_SIDE_FACE
						   gap = parZDash-node2Z - parDia*0.5;
						   if(gap < 0)
						   {
			  				    if (gap < -gapLimit)
							    {
									errorFile<<" DISK_SIDE_FACE "<<gap/id->lengthFactor<<endl;
									cout<<" DISK_SIDE_FACE "<<gap/id->lengthFactor<<endl;
								}
								//errorFile<<" DISK_SIDE_FACE 2 Overlap  "<<endl;
                                
                                unitVec[0] = 0.0;
								unitVec[1] = 0.0;
								unitVec[2] = 1.0;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								//unit vectors for the hole contact point
								//contactVector[0] = -radius*(*unitVec);
								//contactVector[1] = -radius*(*(unitVec+1));
								//contactVector[2] = -radius*(*(unitVec+2));
								jpRVec[0] = parXDash;
								jpRVec[1] = parYDash;
								jpRVec[2] = 0.0;
								
								//errorFile<<" U Vec  "<<unitVec[0]<<"  "<<unitVec[1]<<"  "<<unitVec[2]<<endl;
								isContacted = 1;
								boundary = -8;
								torqueRadius = pXY; //distance to the contact point
								compType = 1; //for torque calculation  
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
								if (SHAFT[parIndex])
								{
									checkBoundary(-2);
								}
								return isContacted;
						   }
						   checkBoundary(-8); //check for disk side face
						   break;
						   
						   case 3: //DISK_EDGE
						   jpX = parXDash*radius/pXY;
						   jpY = parYDash*radius/pXY;
						   jpZ = node2Z;

						   gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar;
						   if(gap < 0)
						   {
				  				if (gap < -gapLimit)
							    {
									errorFile<<" DISK_EDGE "<<gap/id->lengthFactor<<endl;
									cout<<" DISK_EDGE "<<gap/id->lengthFactor<<endl;
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
								
								boundary = -5;
								torqueRadius = 0.5*componentDia[compI]; //distance to the contact point
								compType = 1; //for torque calculation  
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);

								isContacted = 1;
								return isContacted;
						   }
						   checkBoundary(-5); //check for disk edge
						   break;
						   
						   case 4: //DISK_EDGE
						   jpX = parXDash*radius/pXY;
						   jpY = parYDash*radius/pXY;
						   jpZ = node1Z;

						   gap = sqrt(pow(parXDash-jpX,2)+pow(parYDash-jpY,2)+pow(parZDash-jpZ,2)) - rStar;
						   if(gap < 0)
						   {
				  				if (gap < -gapLimit)
							    {
									errorFile<<" DISK_EDGE "<<gap/id->lengthFactor<<endl;
									cout<<" DISK_EDGE "<<gap/id->lengthFactor<<endl;
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
								//unit vectors for the hole contact point
								//contactVector[0] = -radius*(*unitVec);
								//contactVector[1] = -radius*(*(unitVec+1));
								//contactVector[2] = -radius*(*(unitVec+2));
								jpRVec[0] = jpX;     
								jpRVec[1] = jpY; 
								jpRVec[2] = 0.0; 
								
								boundary = -5;
								torqueRadius = 0.5*componentDia[compI]; //distance to the contact point
								compType = 1; //for torque calculation  
	
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
							
								isContacted = 1;
								return isContacted;
						   }
						   checkBoundary(-5); //check for disk edge
						   break;
						   
						   case 5: //DISK_OUTERFACE
						   gap = pXY - radius - rStar;
						   if(gap < 0)
						   {
				  				if (gap < -gapLimit)
							    {
									errorFile<<" DISK_OUTERFACE "<<gap/id->lengthFactor<<endl;
									cout<<" DISK_OUTERFACE "<<gap/id->lengthFactor<<endl;
								}
								unitVec[0] = parXDash/pXY;
								unitVec[1] = parYDash/pXY;
								unitVec[2] = 0.0;
								//unit vectors for the particle
								ipRVec[0] = -0.5*parDia*(*unitVec);
								ipRVec[1] = -0.5*parDia*(*(unitVec+1));
								ipRVec[2] = -0.5*parDia*(*(unitVec+2));
								//cout<< "ipRVec[1]  "<<ipRVec[1]<<endl;
								//unit vectors for the hole contact point
								//contactVector[0] = -radius*(*unitVec);
								//contactVector[1] = -radius*(*(unitVec+1));
								//contactVector[2] = -radius*(*(unitVec+2));
								jpRVec[0] = radius*(*unitVec);     
								jpRVec[1] = radius*(*(unitVec+1)); 
								jpRVec[2] = radius*(*(unitVec+2)); 
								
								//errorFile<<" U VEC "<<unitVec[0]<<"   "<<unitVec[1]<<endl;
								boundary = -4;
								torqueRadius = 0.5*componentDia[compI]; //distance to the contact point
								compType = 1; //for torque calculation  
	
								shaftAngVel[2] = componentAngVelocity[compI];
								findContactForce(1);
								isContacted = 1;
                                return isContacted;
								//rotatingSpeed = rotSpeed;
						    }
						    checkBoundary(-4); //check for disk outer face
						   //cout<<" gap  "<<gap<<endl;
						  
					 break;
			}//end of switch
		}//isInside
	}//particle close enough
	//checkAllBoundary();
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
			errorFile<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" SHAFT"<<endl;
			cout<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" SHAFT"<<endl;
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
		compType = 2; //for torque calculation  
	
		shaftAngVel[2] = componentAngVelocity[compI];
	    findContactForce(1);
		if (DRUM_SIDE_FACE[parIndex])
		{
			checkBoundary(-3);
		}

		if (DISK_SIDE_FACE[parIndex])
		{
			checkBoundary(-8);
		}
		isContacted = 1;
		return isContacted;
	}
	checkBoundary(-2); //check for shaft itself
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
			errorFile<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" DRUM"<<endl;
			cout<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" DRUM"<<endl;
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
		boundary = -1;
		shaftAngVel[2] = componentAngVelocity[compI];
		findContactForce(1);
        //errorFile<<" DRUM "<<shaftAngVel[2]<<endl;
		if (DRUM_SIDE_FACE[parIndex])
		{
			checkBoundary(-3); //check for end wall
		}

		isContacted = 1;
		return isContacted;
	}
	checkBoundary(-1); //check for drum
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
				errorFile<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" DRUM SIDEFACE"<<endl;
				cout<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" DRUM SIDEFACE"<<endl;
			}
			//contactString = "CONTACT_WITH_DRUM_SIDEFACE";
			ipRVec[0] = -0.5*parDia*(*unitVec);
			ipRVec[1] = -0.5*parDia*(*(unitVec+1));
			ipRVec[2] = -0.5*parDia*(*(unitVec+2));
			jpRVec[0] = parX;
			jpRVec[1] = parY;
			jpRVec[2] = 0.0;
			boundary = -3;
			shaftAngVel[2] = componentAngVelocity[compI];
			findContactForce(1);
			if (DRUM_OUTER_FACE[parIndex])
			{
				checkBoundary(-1);
			}
			if (SHAFT[parIndex])
			{
				checkBoundary(-2);
			}
			isContacted = 1;
			return isContacted;
		}
	}
	checkBoundary(-3);
	return isContacted;
};

int Mill::findSecondEndWallContact(int compI)
{
	isContacted = 0;
	sfCoff = id->sfricpWall; 
	//sfCoff = 0.0; //no sliding friction between particle and second end wall

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
				errorFile<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" DRUM SIDEFACE"<<endl;
				cout<<" par "<<particleNo[parIndex]<<"  "<<gap/id->lengthFactor<<" DRUM SIDEFACE"<<endl;
			}
			//contactString = "CONTACT_WITH_DRUM_SIDEFACE";
				
			ipRVec[0] = -0.5*parDia*(*unitVec);
			ipRVec[1] = -0.5*parDia*(*(unitVec+1));
			ipRVec[2] = -0.5*parDia*(*(unitVec+2));
			jpRVec[0] = parX;
			jpRVec[1] = parY;
			jpRVec[2] = 0.0;
			boundary = -3;
			shaftAngVel[2] = componentAngVelocity[compI];
			findContactForce(1);
			if (DRUM_OUTER_FACE[parIndex])
			{
				checkBoundary(-1);
			}
			if (SHAFT[parIndex])
			{
				checkBoundary(-2);
			}
			isContacted = 1;
			return isContacted;
		}
	}
	checkBoundary(-3);
	return isContacted;
};

bool Mill::checkData(int showKE, int parI, int arrI, double x,double y,double z,double cF, 
				double iE,double sE, double gap, double kE,double tkE, double tW)
{
	bool found = false; // To check contact is existing or not
	//--- Finding cell index for each collision point -------------
	x = x + 0.5*id->chamberInnerDia;
	y = y + 0.5*id->chamberInnerDia;
	int ii = (int)floor(x/dx);
	int jj = (int)floor(y/dy);
	int kk = (int)floor(z/dz);
	if (ii < 0)
	{
		ii =0;
	}
	else if(ii > iIndex-1)
	{
		ii = iIndex-1;
	}
        if (jj < 0)
        {
                jj =0;
        }
        else if(jj > jIndex-1)
        {
                jj = jIndex-1;
        }
        if(kk < 0)
	{
		kk = 0;
	}
	else if(kk > kIndex-1)
        {
                kk = kIndex-1;
        }
	int cellIndex = ii + (jj*iIndex) + (kk*iIndex*jIndex);
	if(cellIndex > iIndex*jIndex*kIndex)
	{
		cout<<ii<<"   "<<jj<<"   "<<kk<<endl;
		exit(0);
	}
	//--------------------------------------------------------------
	for (int i=0; i<noOfContacts[parI] ; i++)
	{
		if (contactList[parI][i] == arrI) // Contact  exists
		{	
			found = true;
			shearEnergyDissipation[parI][i] = shearEnergyDissipation[parI][i] + sE;
			cellShearEnergy[cellIndex]  = cellShearEnergy[cellIndex] + sE;
			if (-gap > previousOverlap[parI][i])
			{
				impactEnergyDissipation[parI][i] = impactEnergyDissipation[parI][i] + iE;
				cellImpactEnergy[cellIndex]  = cellImpactEnergy[cellIndex] + iE;
				previousOverlap[parI][i] = -gap;
			}
		}
	}
		
	if (!found) // A new contact
	{ 
		switch(arrI)
		{
			
			case -1:
				DRUM_OUTER_FACE[parI] = true;
				break;
			case -2:
				SHAFT[parI] = true;
				break;
			case -3:
				DRUM_SIDE_FACE[parI] = true;
				break;
			case -8:
				DISK_SIDE_FACE[parI] = true;
				break;
		}

		totalWear[parI] = tW;
		contactList[parI][noOfContacts[parI]] = arrI;
		collisionEnergy[parI][noOfContacts[parI]] = kE;
		tangCollisoinEnergy[parI][noOfContacts[parI]] = tkE;
		cellCollisionEnergy[cellIndex] = cellCollisionEnergy[cellIndex] + kE;
		cellTangCollisionEnergy[cellIndex] = cellTangCollisionEnergy[cellIndex] + tkE;
		cellNoOfContacts[cellIndex] = cellNoOfContacts[cellIndex] + 1;
		noOfContacts[parI] = noOfContacts[parI] + 1;
		totalCollisions[parI] = totalCollisions[parI] + 1;//keeps a track of total number of collisions, particle having up until any time
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
			double sE = shearEnergyDissipation[parI][i];
			double iE = impactEnergyDissipation[parI][i];
			double kE = collisionEnergy[parI][i];
			double tkE = tangCollisoinEnergy[parI][i];
			forceFile.write(reinterpret_cast <const char*>(&iE),sizeof(double));
			forceFile.write(reinterpret_cast <const char*>(&sE),sizeof(double));
			forceFile.write(reinterpret_cast <const char*>(&kE),sizeof(double));
			forceFile.write(reinterpret_cast <const char*>(&tkE),sizeof(double));

			impactEnergyDissipation[parI][i] = 0.0;
			shearEnergyDissipation[parI][i] = 0.0;
			collisionEnergy[parI][i] = 0.0;
			tangCollisoinEnergy[parI][i] = 0.0;
			previousOverlap[parI][i] = 0.0;
			deleted = true;
			break;
			
		}
	}
	
	switch(arrIndex)
	{
			case -1:
				DRUM_OUTER_FACE[parI] = false;
				break;
			case -2:
				SHAFT[parI] = false;
				break;
			case -3:
				DRUM_SIDE_FACE[parI] = false;
				break;
			case -8:
				DISK_SIDE_FACE[parI] = false;
				break;
	}
	return deleted;
};
/***************************************************************
 * This method checks the contact for all the boundaries except
 * particle-particle contacts
 ***************************************************************/
void Mill::checkAllBoundary()
{
	for (int i=0; i<noOfContacts[parIndex] ; i++ )
	{
		if (contactList[parIndex][i] < 0)
		{
			checkBoundary(contactList[parIndex][i]);
		}
	}
};

void Mill::checkBoundary(int prevB)
{
	ipDia = particleDiameter[parIndex];
	rStar = 0.5*particleDiameter[parIndex];
    pXY = sqrt(pow(parX,2)+pow(parY,2));
	
	switch(prevB)
	{
		case -1: // drum
			gap = 0.5*id->chamberInnerDia - pXY - rStar;

			if (gap > collisionGapCoeff)
			{
				if (deleteContact(parIndex, prevB))
				{
					//boundary = prevB;
					//writeCollisionFile();
				};
			}
			break;

		case -2: // Shaft
			gap = pXY - 0.5*shaftDia - rStar;
			if (gap > collisionGapCoeff)
			{
				if (deleteContact(parIndex, prevB))
				{
					//boundary = prevB;
					//writeCollisionFile();
				};
			}
			break;

		case -3: // drum side face
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
				if (deleteContact(parIndex, prevB))
				{
					//boundary = prevB;
					//writeCollisionFile();
				};
			}
			break;
		
		case -4: // disk outer face
			gap = pXY - 0.5*diskDia - rStar;
			if (gap > collisionGapCoeff)
			{
				if (deleteContact(parIndex, prevB))
				{
					//boundary = prevB;
					//writeCollisionFile();
				};
			}
			break;

		case -5: // disk outer line
			jpX = (parX*diskDia/2.0)/pXY;//X coordinate on disk outline//
		
			jpY = (parY*diskDia/2.0)/pXY;//Y coordinate on disk outline//
			if(parZ - nearestDiskZ > 0)
				jpZ = diskThick/2.0;//Z coordinate on disk outline//
			else
				jpZ = -diskThick/2.0;//Z coordinate on disk outline//
			ijVec[0] = parX-jpX;
			ijVec[1] = parY-jpY;
			ijVec[2] = parZ - nearestDiskZ - jpZ;
			unitVector(ijVec);
			unitVec[0] = r1;
			unitVec[1] = r2;
			unitVec[2] = r3;
			
			gap = sqrt(pow(parX-jpX,2)+pow(parY-jpY,2)+pow(parZ - nearestDiskZ - jpZ,2)) - rStar;
			if (gap > collisionGapCoeff)
			{
				if (deleteContact(parIndex, prevB))
				{
					//boundary = prevB;
					//writeCollisionFile();
				};
			}
			break;

		case -6: // hole face
			xx = parX - nearestHoleX;
			yy = parY - nearestHoleY;
			zz = 0.0;
			xy = sqrt(pow(xx,2)+pow(yy,2));
			gap = id->holeDia*0.5 - xy - rStar;
					//errorFile<<" HOLE "<<nearestHoleX<<"  "<<nearestHoleY<<endl;
			if (gap > collisionGapCoeff)
			{
				if (deleteContact(parIndex, prevB))
				{
                    //writeCollisionFile();
				};
			}
			break;

		case -7: // hole line
			xx = parX - nearestHoleX;  
			yy = parY - nearestHoleY;  
			zz = 0.0;
			xy = sqrt(pow(xx,2)+pow(yy,2));
			unitVec[0] = xx/xy;//always towards the centre from the contact point//
			unitVec[1] = yy/xy;
			unitVec[2] = 0.0;
			jpX = nearestHoleX +id->holeDia/2.0*unitVec[0];
			jpY = nearestHoleY +id->holeDia/2.0*unitVec[1];
			if(parZ-nearestDiskZ > 0)
				jpZ = diskThick/2.0;//Z coordinate on disk outline//
			else
				jpZ = -diskThick/2.0;//Z coordinate on disk outline//
			
			ijVec[0] = parX-jpX;
			ijVec[1] = parY-jpY;
			ijVec[2] = parZ-nearestDiskZ-jpZ;
			//cout<<" ijVec "<<ijVec[0]<<"  "<<ijVec[1]<<"  "<<ijVec[2];
			unitVector(ijVec);
			unitVec[0] = r1;  
			unitVec[1] = r2;  
			unitVec[2] = r3;  
			temp = sqrt(pow(parX-jpX,2)+pow(parY-jpY,2)+pow(parZ-nearestDiskZ-jpZ,2));

			gap =  temp - rStar;
			if (gap > collisionGapCoeff)
			{
				if (deleteContact(parIndex, prevB))
				{
					//boundary = prevB;
					//writeCollisionFile();
				};
			}
			break;

		case -8: // disk side face
			unitVec[0] = 0.0;
			unitVec[1] = 0.0;
			if(parZ > nearestDiskZ)
			{
				unitVec[2] =1.0;
			}	
			else
			{
				//ijVec.push_back(-1);
				unitVec[2] = -1.0;
			}
				
			gap = fabs(parZ - nearestDiskZ) - (diskThick/2.0+parDia/2.0);
			if (gap > collisionGapCoeff)
			{
				if (deleteContact(parIndex, prevB))
				{
					//boundary = prevB;
					//writeCollisionFile();
				};
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
		double velArray[3];
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
		projVec(cntPntVel,unitVec,0); //project relative velocity vector
		tangVel = sqrt(r1*r1+r2*r2+r3*r3);

		normDampC = id->pdNormDampC;
		calculateForce();
		//---Torque calculation is only for the components 
		//--- which are mounted on the shaft 
		if(nrmForce < 0)
        {
            nrmForce = 0.0;
        }
        double forceX = -unitVec[0]*nrmForce;
		double forceY = -unitVec[1]*nrmForce;
		double forceZ = -unitVec[2]*nrmForce;

        double resultantF = sqrt(pow(forceX,2)+pow(forceY,2));
        switch(boundary)
        {
            
            case -2: //shaft
            //errorFile<<" shaft "<<endl;
            torque = torque + resultantF*torqueRadius*sfCoff;
            break;
            
            case -4: //disk outer face
            //errorFile<<" disk outer fac "<<endl;
            torque = torque + resultantF*torqueRadius*sfCoff;
            break;
            
            case -5: //disk edge
            //errorFile<<" disk edge "<<endl;
            torque = torque + resultantF*torqueRadius*sfCoff;
            break;
            
            case -6: //hole face
            //errorFile<<" hole face "<<endl;
            torque = torque - forceX*collPointY + forceY*collPointX;
            //torque = torque + resultantF*torqueRadius*sfCoff;
            break;
            
            case -7: //hole edge
            //errorFile<<" hole edge "<<endl;
            torque = torque - forceX*collPointY + forceY*collPointX;
            //torque = torque + resultantF*torqueRadius*sfCoff;
            break;
            
            case -8: //disk side face
            //errorFile<<" disk side fac "<<endl;
            torque = torque + fabs(forceZ)*torqueRadius*sfCoff;
            //torque = torque + resultantF*torqueRadius*sfCoff;
            break;
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
