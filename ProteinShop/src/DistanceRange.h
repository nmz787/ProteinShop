/***************************************************************************************
*cr								
*cr					Copyright (c) 2004, The Regents of the 
*cr	University of California, through Lawrence Berkeley National Laboratory, 
*cr	Univ. of Calif. at Davis, and Lawrence Livermore National Laboratory 
*cr	(subject to receipt of any required approvals from U.S. Dept. of Energy).  
*cr							All rights reserved.
*cr		
*cr		Please see the accompanying LICENSE file for further information.
*cr
***************************************************************************************/

/***********************************************************************
DistanceRange - Class to store ranges of distances.
***********************************************************************/

#ifndef DISTANCERANGE_INCLUDED

#include "MDGeometry.h"

namespace MD {

class DistanceRange
{
    /* Elements: */
    double min,max; // Original distance range
    double min2,max2; // Squared distance range
    
    /* Constructors and destructors: */
    public:
    DistanceRange(double sMin,double sMax)
        :min(sMin),max(sMax),min2(min*min),max2(max*max)
    {
    };
    
    /* Methods: */
    double getMin(void) const
    {
        return min;
    };
    double getMax(void) const
    {
        return max;
    };
    bool isInRange(double distance) const // Checks if a distance is inside the range
    {
        return min<=distance&&distance<=max;
    };
    bool isSqrInRange(double distance2) const // Checks if a squared distance is inside the range
    {
        return min2<=distance2&&distance2<=max2;
    };
    bool areInRange(const Point& p1,const Point& p2) const // Checks if two points are inside distance range
    {
        double distance2=double(Geometry::sqrDist(p1,p2));
        return min2<=distance2&&distance2<=max2;
    };
};

}

#endif
