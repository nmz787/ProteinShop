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
MDGeometry - Basic definitions for geometry in molecular dynamics.
***********************************************************************/

#ifndef MDGEOMETRY_INCLUDED
#define MDGEOMETRY_INCLUDED

#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/Ray.h>

namespace MD {

typedef double Scalar; // Basic scalar type for geometry representations
typedef Geometry::Point<Scalar,3> Position; // Type for atom positions
typedef Geometry::Point<Scalar,3> Point; // Type for generic points
typedef Geometry::Vector<Scalar,3> Vector; // Type for generic vectors
typedef Geometry::Ray<Scalar,3> Ray; // Type for rays (used for picking and such)
}

#endif
