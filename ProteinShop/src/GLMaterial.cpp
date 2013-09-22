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

#include "GLMaterial.h"

/***************************
Methods of class GLMaterial:
***************************/

GLMaterial::GLMaterial(void)
    :ambient(0.2f,0.2f,0.2f,1.0f),
     diffuse(0.8f,0.8f,0.8f,1.0f),
     specular(0.0f,0.0f,0.0f,1.0f),
     shininess(0.0f),
     emission(0.0f,0.0f,0.0f,1.0f)
{
}

GLMaterial::GLMaterial(const GLMaterial::Color& sAmbientDiffuse)
    :ambient(sAmbientDiffuse),
     diffuse(sAmbientDiffuse),
     specular(0.0f,0.0f,0.0f,1.0f),
     shininess(0.0f),
     emission(0.0f,0.0f,0.0f,1.0f)
{
}

GLMaterial::GLMaterial(const GLMaterial::Color& sAmbientDiffuse,const GLMaterial::Color& sSpecular,GLfloat sShininess)
    :ambient(sAmbientDiffuse),
     diffuse(sAmbientDiffuse),
     specular(sSpecular),
     shininess(sShininess),
     emission(0.0f,0.0f,0.0f,1.0f)
{
}

GLMaterial::GLMaterial(const GLMaterial::Color& sAmbient,const GLMaterial::Color& sDiffuse,const GLMaterial::Color& sSpecular,GLfloat sShininess)
    :ambient(sAmbient),
     diffuse(sDiffuse),
     specular(sSpecular),
     shininess(sShininess),
     emission(0.0f,0.0f,0.0f,1.0f)
{
}

GLMaterial::GLMaterial(const GLMaterial::Color& sAmbientDiffuse,const GLMaterial::Color& sSpecular,GLfloat sShininess,const GLMaterial::Color& sEmission)
    :ambient(sAmbientDiffuse),
     diffuse(sAmbientDiffuse),
     specular(sSpecular),
     shininess(sShininess),
     emission(sEmission)
{
}

GLMaterial::GLMaterial(const GLMaterial::Color& sAmbient,const GLMaterial::Color& sDiffuse,const GLMaterial::Color& sSpecular,GLfloat sShininess,const GLMaterial::Color& sEmission)
    :ambient(sAmbient),
     diffuse(sDiffuse),
     specular(sSpecular),
     shininess(sShininess),
     emission(sEmission)
{
}

void glMaterial(GLenum face,const GLMaterial& material)
{
    glMaterial(face,GL_AMBIENT,material.ambient);
    glMaterial(face,GL_DIFFUSE,material.diffuse);
    glMaterial(face,GL_SPECULAR,material.specular);
    glMaterialf(face,GL_SHININESS,material.shininess);
    glMaterial(face,GL_EMISSION,material.emission);
}
