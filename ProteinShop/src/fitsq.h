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
FITSQ - Fortran procedure to align corresponding point sets.  The
correspondences are taken from the order of the coordinate arrays.  To
link, be sure to put the g2c library first in the list.
***********************************************************************/


#ifndef FITSQ_INCLUDED
#define FITSQ_INCLUDED


// typedefs from f2c.h (maybe not available on target system)
typedef long int integer;
typedef double doublereal;


extern "C" {

/*
PROCEDURE:  fitsq_
IN:
    rms : one element, value not used
    x   : 3 * nn array of coordinates of first structure to be superimposed
    y   : 3 * nn array of coordinates of second structure to be superimposed
    nn  : number of points in x and y arrays; 1/3 of the array length
          note that the procedure comments are wrong about this value
    t   : 3 elements, value not used
    b   : 9 elements, value not used

OUT:
    rms : root mean square distance
    x   : unchanged
    y   : unchanged
    nn  : unchanged
    t   : 1 * 3 translation vector mapping second point set into the frame of the first
    b   : 3 * 3 rotation matrix mapping second point set into the frame of the first

    Returns 0.
*/
int fitsq_ (doublereal *rms, doublereal *x, doublereal *y,
            integer *nn, doublereal *t, doublereal *b);

};


#endif

