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
GLMaterial - Class to encapsulate OpenGL material properties.
***********************************************************************/

#ifndef GLMATERIAL_INCLUDED
#define GLMATERIAL_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLTypes.h>

class GLMaterial
{
    /* Embedded classes: */
    public:
    typedef GLColor<GLfloat,4> Color; // Type for colors used in materials
    
    /* Elements: */
    private:
    Color ambient; // Ambient color component
    Color diffuse; // Diffuse color component
    Color specular; // Specular color component
    GLfloat shininess; // Specular lighting exponent
    Color emission; // Emissive color component
    
    /* Constructors and destructors: */
    public:
    GLMaterial(void); // Constructs default material
    GLMaterial(const Color& sAmbientDiffuse); // Constructs diffuse material
    GLMaterial(const Color& sAmbientDiffuse,const Color& sSpecular,GLfloat sShininess); // Constructs specular material
    GLMaterial(const Color& sAmbient,const Color& sDiffuse,const Color& sSpecular,GLfloat sShininess); // Constructs specular material with separate ambient color
    GLMaterial(const Color& sAmbientDiffuse,const Color& sSpecular,GLfloat sShininess,const Color& sEmission); // Constructs specular and emissive material
    GLMaterial(const Color& sAmbient,const Color& sDiffuse,const Color& sSpecular,GLfloat sShininess,const Color& sEmission); // Full initialization
    
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
    GLfloat getShininess(void) const
    {
        return shininess;
    };
    void setShininess(GLfloat newShininess)
    {
        shininess=newShininess;
    };
    Color getEmission(void) const
    {
        return emission;
    };
    void setEmission(const Color& newEmission)
    {
        emission=newEmission;
    };
    friend void glMaterial(GLenum face,const GLMaterial& material); // Sets the material for front- and/or backfaces
};

#endif
