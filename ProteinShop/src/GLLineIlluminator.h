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

#ifndef GLLINEILLUMINATOR_INCLUDED
#define GLLINEILLUMINATOR_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLContextData.h>
#include <GLTypes.h>
#include <GLMaterial.h>

class GLLineIlluminator
{
    /* Embedded classes: */
    public:
    typedef GLVector<GLfloat,3> Vector;
    typedef GLColor<GLfloat,4> Color;
    
    private:
    enum MaterialType
    {
        NONE,INTENSITY,RGBA
    };
    
    struct DataItem:public GLContextData::DataItem
    {
        /* Elements: */
        public:
        MaterialType materialType; // Type of current material
        GLuint materialTextureId; // ID of material texture object
        
        /* Constructors and destructors: */
        DataItem(void)
            :materialType(NONE)
        {
            glGenTextures(1,&materialTextureId);
        };
        ~DataItem(void)
        {
            glDeleteTextures(1,&materialTextureId);
        };
    };
    
    /* Elements: */
    private:
    Vector sceneCenter; // Central point used to convert origin points to directions
    bool autoViewDirection; // If true, calculate view direction from GL matrices
    Vector viewDirection; // Caller-provided view direction
    bool autoLightDirection; // If true, calculate light direction from GL light source
    GLenum autoLightIndex; // Index of GL light source used for illumination
    Vector lightDirection; // Caller-provided light direction
    
    /* Constructors and destructors: */
    public:
    GLLineIlluminator(void);
    
    /* Methods: */
    void init(GLContextData& contextData) const; // Initializes the illuminator (once per GL context)
    
    /* Methods to set the current line material: */
    void setMaterial(GLContextData& contextData,GLfloat ambient,GLfloat diffuse,GLfloat specular,GLfloat shininess) const;
    void setMaterial(GLContextData& contextData,const Color& ambient,const Color& diffuse,const Color& specular,GLfloat shininess) const;
    void setMaterial(GLContextData& contextData,const GLMaterial& material) const;
    
    void setSceneCenter(const Vector& newSceneCenter) // Sets the scene center used to convert points to vectors (in current model coordinates)
    {
        sceneCenter=newSceneCenter;
    };
    void enableAutoView(void) // Enables automatic calculation of view direction
    {
        autoViewDirection=true;
    };
    void disableAutoView(void) // Disables automatic calculation of view direction
    {
        autoViewDirection=false;
    };
    void setViewDirection(const Vector& newViewDirection); // Sets a view direction (in current model coordinates) and disables auto view direction
    void enableAutoLight(GLenum lightIndex) // Enables automatic calculation of light direction
    {
        autoLightDirection=true;
        autoLightIndex=lightIndex;
    };
    void disableAutoLight(void) // Disables automatic calculation of light direction
    {
        autoLightDirection=false;
    };
    void setLightDirection(const Vector& newLightDirection); // Sets a light direction (in current model coordinates) and disables auto light direction
    void enableLighting(GLContextData& contextData) const; // Sets up GL state to render illuminated lines
    void disableLighting(GLContextData& contextData) const; // Turns off illuminated line rendering
};

#endif
