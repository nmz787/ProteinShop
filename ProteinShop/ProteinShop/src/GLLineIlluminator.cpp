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
GLLineIlluminator - Class to render illuminated lines.
***********************************************************************/

#undef NONSTANDARD_TEMPLATES
#define NONSTANDARD_TEMPLATES 1

#include <Math/Math.h>
#include <Geometry/ComponentArray.h>
#include <Geometry/Matrix.h>
#include <GLColorOps.h>

#include "GLLineIlluminator.h"

/**********************************
Methods of class GLLineIlluminator:
**********************************/

GLLineIlluminator::GLLineIlluminator(void)
    :sceneCenter(0.0,0.0,0.0),
     autoViewDirection(true),
     autoLightDirection(true),autoLightIndex(GL_LIGHT0)
{
}

void GLLineIlluminator::init(GLContextData& contextData) const
{
    DataItem* dataItem=new DataItem;
    contextData.addDataItem(this,dataItem);
}

void GLLineIlluminator::setMaterial(GLContextData& contextData,GLfloat ambient,GLfloat diffuse,GLfloat specular,GLfloat shininess) const
{
    /* Create a 2D texture map encoding Phong's lighting model: */
    static GLfloat texture[32][32];
    for(int x=0;x<32;++x)
    {
        GLfloat s=2.0f*(GLfloat(x)+0.5f)/32.0f-1.0f;
        GLfloat oneMinusS2=1.0f-s*s;
        GLfloat ambientDiffuse=diffuse;
        ambientDiffuse*=GLfloat(pow(Math::sqrt(oneMinusS2),2.0f));
        ambientDiffuse+=ambient;
        for(int y=0;y<32;++y)
        {
            GLfloat t=2.0f*(GLfloat(y)+0.5f)/32.0f-1.0f;
            GLfloat oneMinusT2=1.0f-t*t;
            GLfloat color=specular;
            color*=GLfloat(pow(Math::abs(Math::sqrt(oneMinusS2*oneMinusT2)-s*t),shininess));
            color+=ambientDiffuse;
            texture[y][x]=color;
        }
    }
    
    /* Get the context data item and upload the created texture: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    dataItem->materialType=INTENSITY;
    glBindTexture(GL_TEXTURE_2D,dataItem->materialTextureId);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glPixelStorei(GL_UNPACK_ROW_LENGTH,0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexImage2D(GL_TEXTURE_2D,0,GL_INTENSITY,32,32,0,GL_LUMINANCE,GL_FLOAT,texture);
    glBindTexture(GL_TEXTURE_2D,0);
}

void GLLineIlluminator::setMaterial(GLContextData& contextData,const GLLineIlluminator::Color& ambient,const GLLineIlluminator::Color& diffuse,const GLLineIlluminator::Color& specular,GLfloat shininess) const
{
    /* Create a 2D texture map encoding Phong's lighting model: */
    static Color texture[32][32];
    for(int x=0;x<32;++x)
    {
        GLfloat s=2.0f*(GLfloat(x)+0.5f)/32.0f-1.0f;
        GLfloat oneMinusS2=1.0f-s*s;
        Color ambientDiffuse=diffuse;
        ambientDiffuse*=GLfloat(pow(Math::sqrt(oneMinusS2),2.0f));
        ambientDiffuse+=ambient;
        for(int y=0;y<32;++y)
        {
            GLfloat t=2.0f*(GLfloat(y)+0.5f)/32.0f-1.0f;
            GLfloat oneMinusT2=1.0f-t*t;
            Color color=specular;
            color*=GLfloat(pow(Math::abs(Math::sqrt(oneMinusS2*oneMinusT2)-s*t),shininess));
            color+=ambientDiffuse;
            texture[y][x]=color;
        }
    }
    
    /* Get the context data item and upload the created texture: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    dataItem->materialType=RGBA;
    glBindTexture(GL_TEXTURE_2D,dataItem->materialTextureId);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glPixelStorei(GL_UNPACK_ROW_LENGTH,0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,32,32,0,GL_RGBA,GL_FLOAT,texture);
    glBindTexture(GL_TEXTURE_2D,0);
}

void GLLineIlluminator::setMaterial(GLContextData& contextData,const GLMaterial& material) const
{
    /* Create a 2D texture map encoding Phong's lighting model: */
    static Color texture[32][32];
    for(int x=0;x<32;++x)
    {
        GLfloat s=2.0f*(GLfloat(x)+0.5f)/32.0f-1.0f;
        GLfloat oneMinusS2=1.0f-s*s;
        Color ambientDiffuse=material.getDiffuse();
        ambientDiffuse*=GLfloat(pow(Math::sqrt(oneMinusS2),2.0f));
        ambientDiffuse+=material.getAmbient();
        for(int y=0;y<32;++y)
        {
            GLfloat t=2.0f*(GLfloat(y)+0.5f)/32.0f-1.0f;
            GLfloat oneMinusT2=1.0f-t*t;
            Color color=material.getSpecular();
            color*=GLfloat(pow(Math::abs(Math::sqrt(oneMinusS2*oneMinusT2)-s*t),material.getShininess()));
            color+=ambientDiffuse;
            texture[y][x]=color;
        }
    }
    
    /* Get the context data item and upload the created texture: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    dataItem->materialType=RGBA;
    glBindTexture(GL_TEXTURE_2D,dataItem->materialTextureId);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glPixelStorei(GL_UNPACK_ROW_LENGTH,0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,32,32,0,GL_RGBA,GL_FLOAT,texture);
    glBindTexture(GL_TEXTURE_2D,0);
}

void GLLineIlluminator::setViewDirection(const GLLineIlluminator::Vector& newViewDirection)
{
    autoViewDirection=false;
    viewDirection=newViewDirection;
    GLfloat viewDirectionLen=Math::sqrt(Math::sqr(viewDirection[0])+Math::sqr(viewDirection[1])+Math::sqr(viewDirection[2]));
    viewDirection[0]/=viewDirectionLen;
    viewDirection[1]/=viewDirectionLen;
    viewDirection[2]/=viewDirectionLen;
}

void GLLineIlluminator::setLightDirection(const GLLineIlluminator::Vector& newLightDirection)
{
    autoLightDirection=false;
    lightDirection=newLightDirection;
    GLfloat lightDirectionLen=Math::sqrt(Math::sqr(lightDirection[0])+Math::sqr(lightDirection[1])+Math::sqr(lightDirection[2]));
    lightDirection[0]/=lightDirectionLen;
    lightDirection[1]/=lightDirectionLen;
    lightDirection[2]/=lightDirectionLen;
}

void GLLineIlluminator::enableLighting(GLContextData& contextData) const
{
    GLenum previousMatrixMode=glGet<GLint>(GL_MATRIX_MODE);
    
    Geometry::Matrix<GLfloat,4,4> modelView;
    if(autoViewDirection||autoLightDirection)
    {
        /* Get the modelview matrix from OpenGL: */
        GLfloat matrixArray[16];
        glGetFloatv(GL_MODELVIEW_MATRIX,matrixArray);
        modelView=Geometry::Matrix<GLfloat,4,4>::fromColumnMajor(matrixArray);
    }
    
    /* Determine the view direction: */
    Geometry::ComponentArray<GLfloat,3> viewDir(viewDirection.getComponents());
    if(autoViewDirection)
    {
        /* Get the projection matrix from OpenGL: */
        GLfloat matrixArray[16];
        glGetFloatv(GL_PROJECTION_MATRIX,matrixArray);
        Geometry::Matrix<GLfloat,4,4> projection=Geometry::Matrix<GLfloat,4,4>::fromColumnMajor(matrixArray);
        
        /* Calculate the view direction from the OpenGL projection and modelview matrices: */
        Geometry::ComponentArray<GLfloat,4> viewPos(0.0f,0.0f,1.0f,0.0f);
        viewPos=viewPos/projection;
        viewPos=viewPos/modelView;
        
        /* Check if it's an orthogonal or perspective projection: */
        if(Math::abs(viewPos[3])<1.0e-8f)
        {
            /* Just copy the view direction: */
            viewDir=viewPos;
        }
        else
        {
            /* Calculate the direction from the view point to the scene center: */
            for(int i=0;i<3;++i)
                viewDir[i]=viewPos[i]/viewPos[3]-sceneCenter[i];
        }
        GLfloat viewDirLen=GLfloat(Geometry::mag(viewDir));
        for(int i=0;i<3;++i)
            viewDir[i]/=viewDirLen;
    }
    
    /* Determine the light direction: */
    Geometry::ComponentArray<GLfloat,3> lightDir(lightDirection.getComponents());
    if(autoLightDirection)
    {
        /* Query the light direction from OpenGL and transform it to model coordinates: */
        Geometry::ComponentArray<GLfloat,4> lightPos;
        glGetLight(autoLightIndex,GL_POSITION,lightPos.getComponents());
        lightPos=lightPos/modelView;
        
        /* Check if it's a directional or point light: */
        if(Math::abs(lightPos[3])<1.0e-8f)
        {
            /* Just copy the light direction: */
            lightDir=lightPos;
        }
        else
        {
            /* Calculate the direction from the light source to the scene center: */
            for(int i=0;i<3;++i)
                lightDir[i]=lightPos[i]/lightPos[3]-sceneCenter[i];
        }
        GLfloat lightDirLen=GLfloat(Geometry::mag(lightDir));
        for(int i=0;i<3;++i)
            lightDir[i]/=lightDirLen;
    }
    
    /* Set up the OpenGL texture matrix: */
    glMatrixMode(GL_TEXTURE);
    glPushMatrix();
    GLfloat matrix[4][4];
    for(int j=0;j<3;++j)
    {
        matrix[j][0]=lightDir[j];
        matrix[j][1]=viewDir[j];
        matrix[j][2]=0.0f;
        matrix[j][3]=0.0f;
    }
    matrix[3][0]=1.0f;
    matrix[3][1]=1.0f;
    matrix[3][2]=0.0f;
    matrix[3][3]=2.0f;
    glLoadMatrixf((const GLfloat*)matrix);
    
    /* Set the OpenGL rendering mode: */
    glPushAttrib(GL_TEXTURE_BIT);
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    glBindTexture(GL_TEXTURE_2D,dataItem->materialTextureId);
    glEnable(GL_TEXTURE_2D);
    if(dataItem->materialType==INTENSITY)
        glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
    else
        glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    
    /* Clean up: */
    glMatrixMode(previousMatrixMode);
}

void GLLineIlluminator::disableLighting(GLContextData& contextData) const
{
    GLenum previousMatrixMode=glGet<GLint>(GL_MATRIX_MODE);
    
    /* Reset the texture matrix: */
    glMatrixMode(GL_TEXTURE);
    glPopMatrix();
    
    /* Reset the OpenGL rendering mode: */
    glPopAttrib();
    
    /* Clean up: */
    glMatrixMode(previousMatrixMode);
}
