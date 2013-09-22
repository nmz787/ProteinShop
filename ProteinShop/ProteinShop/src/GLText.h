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
GLText - Class to encapsulate OpenGL text properties.
***********************************************************************/

#ifndef GLTEXT_H
#define GLTEXT_H

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <FL/Fl.H>
#include <FL/gl.h>

class GLText
{
   /* Elements: */
    private:
	int pw,ph;
   	Fl_Font fontType;
	int fontSize;
	GLfloat fontColor[3];
	GLfloat fontScale[2]; 

	/* Constructors and destructors: */
    public:
   	GLText(void); // Constructs text
	~GLText(void) {}; 

	char fontbuffer[20];

    /* Methods: */
	void setFontType(Fl_Font font);
	void setFontSize(int size);
	void resetFontSize();
	void setFontColor(float r, float g, float b);
	void setFontScale(float x, float y);
	
	/* Fltk font */
	void drawFreeString(const char *string, float x, float y);
	void drawFreeString(const char *string);
	void drawFlatString(const char *string, float x, float y);

	void setWinVars(int w, int h);
};

#endif
