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
GLColorMap - Class for mapping scalar values to RGBA colors.
***********************************************************************/

#include <math.h>
#include <string.h>
#include <FileIO.h>
#include <Endianness.h>

#include "GLColorMap.h"

/***************************
Methods of class GLColorMap:
***************************/

GLColorMap::GLColorMap(GLsizei sNumEntries,GLdouble sMin,GLdouble sMax)
    :numEntries(sNumEntries),entries(new Color[numEntries]),min(sMin),max(sMax),
     factor(GLdouble(numEntries-1)/(max-min)),offset(min*factor)
{
}

GLColorMap::GLColorMap(GLenum type,GLfloat alphaMax,GLfloat alphaGamma,GLdouble sMin,GLdouble sMax)
    :numEntries(256),entries(new Color[numEntries]),min(sMin),max(sMax),
     factor(GLdouble(numEntries-1)/(max-min)),offset(min*factor)
{
    /* Create the palette colors: */
    if(type&GREYSCALE)
    {
        for(int i=0;i<256;++i)
        {
            entries[i][0]=GLfloat(i)/255.0f;
            entries[i][1]=GLfloat(i)/255.0f;
            entries[i][2]=GLfloat(i)/255.0f;
        }
    }
    else if(type&RAINBOW)
    {
        for(int i=0;i<256;++i)
        {
            /* Create rainbow colors: */
            GLdouble rad=GLdouble(i)*(2.0*M_PI/256.0);
            if(rad<=2.0*M_PI/3.0)
                entries[i][0]=cos(0.75*rad);
            else if(rad>=4.0*M_PI/3.0)
                entries[i][0]=cos(0.75*(2.0*M_PI-rad));
            else
                entries[i][0]=0.0;
            entries[i][1]=sin(0.75*rad);
            if(entries[i][1]<0.0)
                entries[i][1]=0.0;
            entries[i][2]=sin(0.75*(rad-2.0*M_PI/3.0));
            if(entries[i][2]<0.0)
                entries[i][2]=0.0;
        }
    }

    /* Create the palette opacities: */
    if(type&CONSTANT_ALPHA)
    {
        for(int i=0;i<256;++i)
            entries[i][3]=alphaMax;
    }
    else if(type&RAMP_ALPHA)
    {
        for(int i=0;i<256;++i)
            entries[i][3]=alphaMax*pow(double(i)/255.0,double(alphaGamma));
    }
}

GLColorMap::GLColorMap(GLsizei sNumEntries,const GLColorMap::Color* sEntries,GLdouble sMin,GLdouble sMax)
    :numEntries(sNumEntries),entries(new Color[numEntries]),min(sMin),max(sMax),
     factor(GLdouble(numEntries-1)/(max-min)),offset(min*factor)
{
    memcpy(entries,sEntries,numEntries*sizeof(Color));
}

GLColorMap::GLColorMap(const char* filename,GLdouble sMin,GLdouble sMax)
    :numEntries(256),entries(new Color[256]),min(sMin),max(sMax),
     factor(GLdouble(numEntries-1)/(max-min)),offset(min*factor)
{
    /* Load the color entries from file: */
    FILE* paletteFile=fopen(filename,"rb");
    if(paletteFile!=0)
    {
        fread(entries,256,paletteFile);
        fclose(paletteFile);
        #if __BIG_ENDIAN!=1234
        for(int i=0;i<256;++i)
            for(int j=0;j<4;++j)
                swapEndianness(entries[i][j]);
        #endif
    }
}

GLColorMap::GLColorMap(const GLColorMap& source)
    :numEntries(source.numEntries),entries(new Color[numEntries]),min(source.min),
     max(source.max),factor(source.factor),offset(source.offset)
{
    memcpy(entries,source.entries,numEntries*sizeof(Color));
}

GLColorMap& GLColorMap::operator=(const GLColorMap& source)
{
    if(this!=&source)
    {
        /* Reallocate the entry array if the number of entries has to change: */
        if(numEntries!=source.numEntries)
        {
            delete[] entries;
            numEntries=source.numEntries;
            entries=new Color[numEntries];
        }

        memcpy(entries,source.entries,numEntries*sizeof(Color));
        min=source.min;
        max=source.max;
        factor=source.factor;
        offset=source.offset;
    }

    return *this;
}

GLColorMap::~GLColorMap(void)
{
    delete[] entries;
}

GLColorMap& GLColorMap::load(const char* filename)
{
    /* Reallocate the entry array if the number of entries has to change: */
    if(numEntries!=256)
    {
        delete[] entries;
        numEntries=256;
        entries=new Color[numEntries];
        factor=GLdouble(numEntries-1)/(max-min);
        offset=min*factor;
    }

    /* Load the color entries from file: */
    FILE* paletteFile=fopen(filename,"rb");
    if(paletteFile!=0)
    {
        fread(entries,256,paletteFile);
        fclose(paletteFile);
        #if __BIG_ENDIAN!=1234
        for(int i=0;i<256;++i)
            for(int j=0;j<4;++j)
                swapEndianness(entries[i][j]);
        #endif
    }

    return *this;
}

GLColorMap& GLColorMap::setColors(int newNumEntries,const GLColorMap::Color* newEntries)
{
    /* Reallocate the entry array if the number of entries has to change: */
    if(numEntries!=newNumEntries)
    {
        delete[] entries;
        numEntries=newNumEntries;
        entries=new Color[numEntries];
        factor=GLdouble(numEntries-1)/(max-min);
        offset=min*factor;
    }

    /* Copy the new entry array: */
    memcpy(entries,newEntries,numEntries*sizeof(Color));

    return *this;
}

void GLColorMap::save(const char* filename) const
{
    /* We only save 256-entry maps! */
    if(numEntries!=256)
        return;
    
    FILE* output=fopen(filename,"wb");
    if(output==0)
        return;
    #if __BIG_ENDIAN!=1234
    GLColor<GLfloat,4> tempEntries[256];
    for(int i=0;i<256;++i)
    {
        tempEntries[i]=entries[i];
        for(int j=0;j<4;++j)
            swapEndianness(tempEntries[i][j]);
    }
    fwrite(tempEntries,numEntries,output);
    #else
    fwrite(entries,numEntries,output);
    #endif
    fclose(output);
}

GLColorMap& GLColorMap::setScalarRange(GLdouble newMin,GLdouble newMax)
{
    min=newMin;
    max=newMax;
    factor=GLdouble(numEntries-1)/(max-min);
    offset=min*factor;

    return *this;
}

GLColorMap& GLColorMap::changeTransparency(GLfloat gamma)
{
    /* Change the transparencies (not opacities!): */
    for(int i=0;i<numEntries;++i)
        entries[i][3]=GLfloat(1.0f-pow(1.0f-entries[i][3],gamma));

    return *this;
}

GLColorMap& GLColorMap::premultiplyAlpha(void)
{
    /* Premultiply each entry: */
    for(int i=0;i<numEntries;++i)
    {
        for(int j=0;j<3;++j)
            entries[i][j]*=entries[i][3];
    }

    return *this;
}

GLColorMap::Color GLColorMap::mapColor(GLdouble scalar) const
{
    /* Check for out-of-bounds arguments: */
    if(scalar<=min)
        return entries[0];
    else if(scalar>=max)
        return entries[numEntries-1];

    /* Calculate the base map index: */
    scalar=scalar*factor-offset;
    GLint index=GLint(floor(scalar));
    if(index==numEntries-1)
        --index;
    scalar-=GLdouble(index);
    Color result;
    result[0]=GLfloat(entries[index][0]*(1.0f-scalar)+entries[index+1][0]*scalar);
    result[1]=GLfloat(entries[index][1]*(1.0f-scalar)+entries[index+1][1]*scalar);
    result[2]=GLfloat(entries[index][2]*(1.0f-scalar)+entries[index+1][2]*scalar);
    result[3]=GLfloat(entries[index][3]*(1.0f-scalar)+entries[index+1][3]*scalar);

    return result;
}
