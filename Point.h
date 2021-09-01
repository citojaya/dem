/**************************************************************
 * Point.h
 * Copyright: Copyright (c) 2003
 * SIMPAS - University of New South Wales, Sydney
 * @author Chandana Jayasundara 		
 * @version 1.0
 * 
 **************************************************************/

#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<iostream>
using namespace std;
#ifndef POINT_H
#define POINT_H

class Point
{
	public:
	Point(double x, double y, double z, double pt, double pt2);
	~Point();
	double *getPoint();
    double pointData,pointData2;
	
	private:
		double *point;
};

#endif
