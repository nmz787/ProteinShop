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
GLLightSource - Class to encapsulate OpenGL light sources.
***********************************************************************/

#include "GLLightSource.h"

/******************************
Methods of class GLLightSource:
******************************/

GLLightSource::GLLightSource(void)
    :ambient(0.0f,0.0f,0.0f,1.0f),
     diffuse(1.0f,1.0f,1.0f,1.0f),
     specular(1.0f,1.0f,1.0f,1.0f),
     position(0.0f,0.0f,1.0f,0.0f),
     spotDirection(0.0f,0.0f,-1.0f),
     spotCutoff(180.0f),
     spotExponent(0.0f),
     constantAttenuation(1.0f),
     linearAttenuation(0.0f),
     quadraticAttenuation(0.0f)
{
}

GLLightSource::GLLightSource(const GLLightSource::Color& sColor,const GLLightSource::Vector& sPosition,GLfloat sConstantAttenuation,GLfloat sLinearAttenuation,GLfloat sQuadraticAttenuation)
    :ambient(0.0f,0.0f,0.0f,1.0f),
     diffuse(sColor),
     specular(sColor),
     position(sPosition),
     spotDirection(0.0f,0.0f,-1.0f),
     spotCutoff(180.0f),
     spotExponent(0.0f),
     constantAttenuation(sConstantAttenuation),
     linearAttenuation(sLinearAttenuation),
     quadraticAttenuation(sQuadraticAttenuation)
{
}

GLLightSource::GLLightSource(const GLLightSource::Color& sColor,const GLLightSource::Vector& sPosition,const GLLightSource::Vector& sSpotDirection,GLfloat sSpotCutoff,GLfloat sSpotExponent,GLfloat sConstantAttenuation,GLfloat sLinearAttenuation,GLfloat sQuadraticAttenuation)
    :ambient(0.0f,0.0f,0.0f,1.0f),
     diffuse(sColor),
     specular(sColor),
     position(sPosition),
     spotDirection(sSpotDirection),
     spotCutoff(sSpotCutoff),
     spotExponent(sSpotExponent),
     constantAttenuation(sConstantAttenuation),
     linearAttenuation(sLinearAttenuation),
     quadraticAttenuation(sQuadraticAttenuation)
{
}

GLLightSource::GLLightSource(const GLLightSource::Color& sAmbient,const GLLightSource::Color& sDiffuse,const GLLightSource::Color& sSpecular,const GLLightSource::Vector& sPosition,const GLLightSource::Vector& sSpotDirection,GLfloat sSpotCutoff,GLfloat sSpotExponent,GLfloat sConstantAttenuation,GLfloat sLinearAttenuation,GLfloat sQuadraticAttenuation)
    :ambient(sAmbient),
     diffuse(sDiffuse),
     specular(sSpecular),
     position(sPosition),
     spotDirection(sSpotDirection),
     spotCutoff(sSpotCutoff),
     spotExponent(sSpotExponent),
     constantAttenuation(sConstantAttenuation),
     linearAttenuation(sLinearAttenuation),
     quadraticAttenuation(sQuadraticAttenuation)
{
}

void glLight(GLenum light,const GLLightSource& lightSource)
{
    glLight(light,GL_AMBIENT,lightSource.ambient);
    glLight(light,GL_DIFFUSE,lightSource.diffuse);
    glLight(light,GL_SPECULAR,lightSource.specular);
    glLight(light,GL_POSITION,lightSource.position);
    glLight(light,GL_SPOT_DIRECTION,lightSource.spotDirection);
    glLightf(light,GL_SPOT_CUTOFF,lightSource.spotCutoff);
    glLightf(light,GL_SPOT_EXPONENT,lightSource.spotExponent);
    glLightf(light,GL_CONSTANT_ATTENUATION,lightSource.constantAttenuation);
    glLightf(light,GL_LINEAR_ATTENUATION,lightSource.linearAttenuation);
    glLightf(light,GL_QUADRATIC_ATTENUATION,lightSource.quadraticAttenuation);
}
