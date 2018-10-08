/**************************************************************
 * Point.cpp
 * This class is used to store finite element data points of 
 * three disks
 * Copyright: Copyright (c) 2003
 * SIMPAS - University of New South Wales, Sydney
 * @author Chandana Jayasundara 		
 * @version 1.0
 * 
 **************************************************************/
#include "Point.h"

Point::Point(double x, double y,double z, double pt, double pt2)
{
	point = new double[3];

	point[0] = x;
	point[1] = y;
	point[2] = z;
    pointData = pt;
    pointData2 = pt2;
};

Point::~Point()
{
	//cout<<"POITN"<<endl;
	//delete []point;
};

double *Point::getPoint()
{
	return point;
};
