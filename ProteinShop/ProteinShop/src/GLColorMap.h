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

#ifndef GLCOLORMAP_INCLUDED
#define GLCOLORMAP_INCLUDED

#include <GLTypes.h>

class GLColorMap
{
    /* Embedded classes: */
    public:
    typedef GLColor<GLfloat,4> Color; // Type of color entries

    enum CreationTypes // Types for automatic palette generation
    {
        GREYSCALE=0x1,RAINBOW=0x2,
        CONSTANT_ALPHA=0x4,RAMP_ALPHA=0x8
    };

    /* Elements: */
    private:
    GLsizei numEntries; // Number of colors in the map
    Color* entries; // Array of RGBA entries
    GLdouble min,max; // The scalar value range
    GLdouble factor,offset; // The scaling factors to map data values to indices

    /* Constructors and destructors: */
    public:
    GLColorMap(void) // Creates an empty color map
        :numEntries(0),entries(0)
    {
    };
    GLColorMap(GLsizei sNumEntries, GLdouble sMin, GLdouble sMax); // Creates an uninitialized color map with the specified number of entries
    GLColorMap(GLenum type,GLfloat alphaMax,GLfloat alphaGamma,GLdouble sMin,GLdouble sMax); // Creates a 256-entry standard color map
    GLColorMap(GLsizei sNumEntries,const Color* sEntries,GLdouble sMin,GLdouble sMax); // Creates a color map from a color array
    GLColorMap(const char* filename,GLdouble sMin,GLdouble sMax); // Loads a 256 entry palette from a file
    GLColorMap(const GLColorMap& source);
    GLColorMap& operator=(const GLColorMap& source);
    ~GLColorMap(void);

    /* Methods: */
    GLColorMap& load(const char* filename); // Loads a 256-entry color map from a file
    GLColorMap& setColors(int newNumEntries,const Color* newEntries); // Sets the color map array directly
    void save(const char* filename) const; // Saves a 256-entry color map to a file
    GLdouble getScalarRangeMin(void) const // Returns minimum of scalar value range
    {
        return min;
    };
    GLdouble getScalarRangeMax(void) const // Returns maximum of scalar value range
    {
        return max;
    };
    GLColorMap& setScalarRange(GLdouble newMin,GLdouble newMax); // Changes the scalar value range
    GLColorMap& changeTransparency(GLfloat gamma); // Applies a gamma function to the transparency values
    GLColorMap& premultiplyAlpha(void); // Converts the colors into premultiplied alpha format for easier compositing
    int getNumEntries(void) const // Returns the number of entries in the map
    {
        return numEntries;
    };
    const Color* getColors(void) const // Returns a pointer to the color entry array
    {
        return entries;
    };
    Color mapColor(GLdouble scalar) const; // Returns the color for a scalar value using linear interpolation
    Color &operator[] (int i) { return entries[i]; }
};

#endif
