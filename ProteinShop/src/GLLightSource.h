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

#ifndef GLLIGHTSOURCE_INCLUDED
#define GLLIGHTSOURCE_INCLUDED

#include <GLTypes.h>

class GLLightSource
{
    /* Embedded classes: */
    public:
    typedef GLVector<GLfloat,4> Vector;
    typedef GLColor<GLfloat,4> Color;
    
    /* Elements: */
    private:
    Color ambient; // Ambient light color
    Color diffuse; // Diffuse light color
    Color specular; // Specular light color
    Vector position; // Point light position or directional light source direction
    Vector spotDirection; // Spot light direction (w component is ignored but should be zero)
    GLfloat spotCutoff; // Spot light cutoff angle
    GLfloat spotExponent; // Spot light attenuation exponent
    GLfloat constantAttenuation; // Point light constant attenuation coefficient
    GLfloat linearAttenuation; // Point light linear attenuation coefficient
    GLfloat quadraticAttenuation; // Point light quadratic attenuation coefficient
    
    /* Constructors and destructors: */
    public:
    GLLightSource(void); // Constructs default light source
    GLLightSource(const Color& sColor,const Vector& sPosition,GLfloat sConstantAttenuation =1.0f,GLfloat sLinearAttenuation =0.0f,GLfloat sQuadraticAttenuation =0.0f); // Sets a monochromatic point or directional light source
    GLLightSource(const Color& sColor,const Vector& sPosition,const Vector& sSpotDirection,GLfloat sSpotCutoff =180.0f,GLfloat sSpotExponent =0.0f,GLfloat sConstantAttenuation =1.0f,GLfloat sLinearAttenuation =0.0f,GLfloat sQuadraticAttenuation =0.0f); // Sets a monochromatic spot light source
    GLLightSource(const Color& sAmbient,const Color& sDiffuse,const Color& sSpecular,const Vector& sPosition,const Vector& sSpotDirection,GLfloat sSpotCutoff,GLfloat sSpotExponent,GLfloat sConstantAttenuation =1.0f,GLfloat sLinearAttenuation =0.0f,GLfloat sQuadraticAttenuation =0.0f); // Full initialization
    
    /* Methods: */
    Color getAmbient(void) const
    {
        return ambient;
    };
    void setAmbient(const Color& newAmbient)
    {
        ambient=newAmbient;
    };
    Color getDiffuse(void) const
    {
        return diffuse;
    };
    void setDiffuse(const Color& newDiffuse)
    {
        diffuse=newDiffuse;
    };
    void setAmbientDiffuse(const Color& newAmbientDiffuse)
    {
        ambient=newAmbientDiffuse;
        diffuse=newAmbientDiffuse;
    };
    Color getSpecular(void) const
    {
        return specular;
    };
    void setSpecular(const Color& newSpecular)
    {
        specular=newSpecular;
    };
    void setColor(const Color& newColor)
    {
        ambient=newColor;
        diffuse=newColor;
        specular=newColor;
    };
    Vector getPosition(void) const
    {
        return position;
    };
    void setPosition(const Vector& newPosition)
    {
        position=newPosition;
    };
    Vector getSpotDirection(void) const
    {
        return spotDirection;
    };
    void setSpotDirection(const Vector& newSpotDirection)
    {
        spotDirection=newSpotDirection;
    };
    GLfloat getSpotCutoff(void) const
    {
        return spotCutoff;
    };
    void setSpotCutoff(GLfloat newSpotCutoff)
    {
        spotCutoff=newSpotCutoff;
    };
    GLfloat getSpotExponent(void) const
    {
        return spotExponent;
    };
    void setSpotExponent(GLfloat newSpotExponent)
    {
        spotExponent=newSpotExponent;
    };
    void setSpotlight(const Vector& newSpotDirection,GLfloat newSpotCutoff,GLfloat newSpotExponent)
    {
        spotDirection=newSpotDirection;
        spotCutoff=newSpotCutoff;
        spotExponent=newSpotExponent;
    };
    GLfloat getConstantAttenuation(void) const
    {
        return constantAttenuation;
    };
    void setConstantAttenuation(GLfloat newConstantAttenuation)
    {
        constantAttenuation=newConstantAttenuation;
    };
    GLfloat getLinearAttenuation(void) const
    {
        return linearAttenuation;
    };
    void setLinearAttenuation(GLfloat newLinearAttenuation)
    {
        linearAttenuation=newLinearAttenuation;
    };
    GLfloat getQuadraticAttenuation(void) const
    {
        return quadraticAttenuation;
    };
    void setQuadraticAttenuation(GLfloat newQuadraticAttenuation)
    {
        quadraticAttenuation=newQuadraticAttenuation;
    };
    void setAttenuation(GLfloat newConstantAttenuation,GLfloat newLinearAttenuation,GLfloat newQuadraticAttenuation)
    {
        constantAttenuation=newConstantAttenuation;
        linearAttenuation=newLinearAttenuation;
        quadraticAttenuation=newQuadraticAttenuation;
    };
    friend void glLight(GLenum light,const GLLightSource& lightSource); // Sets an OpenGL light source
};

#endif
