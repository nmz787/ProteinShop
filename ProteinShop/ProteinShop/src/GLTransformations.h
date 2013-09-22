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
GLTransformations - Wrappers to use Templatized Geometry Library
transformation classes in combination with OpenGL transformation
functions.
***********************************************************************/

#ifndef GLTRANSFORMATIONS_INCLUDED
#define GLTRANSFORMATIONS_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <Geometry/Rotation.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/OrthogonalTransformation.h>
#include <Geometry/AffineTransformation.h>
#include <Geometry/ProjectiveTransformation.h>

template <class ScalarParam>
void glRotate(const Geometry::Rotation<ScalarParam,3>& r);
template <class ScalarParam>
void glLoadMatrix(const Geometry::Rotation<ScalarParam,3>& r);
template <class ScalarParam>
void glMultMatrix(const Geometry::Rotation<ScalarParam,3>& r);
template <class ScalarParam>
void glLoadMatrix(const Geometry::OrthonormalTransformation<ScalarParam,3>& t);
template <class ScalarParam>
void glMultMatrix(const Geometry::OrthonormalTransformation<ScalarParam,3>& t);
template <class ScalarParam>
void glLoadMatrix(const Geometry::OrthogonalTransformation<ScalarParam,3>& t);
template <class ScalarParam>
void glMultMatrix(const Geometry::OrthogonalTransformation<ScalarParam,3>& t);
template <class ScalarParam>
void glLoadMatrix(const Geometry::AffineTransformation<ScalarParam,3>& t);
template <class ScalarParam>
void glMultMatrix(const Geometry::AffineTransformation<ScalarParam,3>& t);
template <class ScalarParam>
void glLoadMatrix(const Geometry::ProjectiveTransformation<ScalarParam,3>& t);
template <class ScalarParam>
void glMultMatrix(const Geometry::ProjectiveTransformation<ScalarParam,3>& t);
template <class ScalarParam>
Geometry::ProjectiveTransformation<ScalarParam,3> glGetMatrix(GLenum whichMatrix);

#if defined(NONSTANDARD_TEMPLATES) && !defined(GLTRANSFORMATIONS_IMPLEMENTATION)
#include <GLTransformations.cpp>
#endif

#endif
