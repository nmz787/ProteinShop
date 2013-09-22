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

#include "GLText.h"
#include <Geometry/Matrix.h>

/***************************
Methods of class GLText:
***************************/
GLText::GLText(void)
:fontType(FL_SCREEN),fontSize(12)
{
	setFontColor(0.0, 0.0, 0.0);
	setFontScale(0.1, 0.1);
}

void GLText::setWinVars(int w, int h)
{
	pw = w;ph = h;
}

void GLText::setFontType(Fl_Font font)
{
	fontType = font;
}

void GLText::setFontSize(int size)
{
	fontSize = size;
	setFontScale(size/100.0, size/100.0);
}

void GLText::resetFontSize()
{
	setFontSize(12);
}

void GLText::setFontColor(float r, float g, float b)
{
	fontColor[0] =r; fontColor[1] =g; fontColor[2] =b;
}

void GLText::setFontScale(float x, float y)
{
	fontScale[0] =x; fontScale[1] =y; 
}

void GLText::drawFlatString(const char *string, float x, float y)
{
	Geometry::Matrix<float,4,4> tempmat;
    glGetFloatv(GL_PROJECTION_MATRIX, (GLfloat *)&tempmat);
    
	glPushMatrix();
   	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
 	glOrtho(0.0f, pw, 0.0f, ph, -1, 1); 
    
	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  	glColor3fv(fontColor);
  	gl_font(fontType, fontSize );
  	gl_draw(string, x, y );

	glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf((GLfloat *)&tempmat);
    glMatrixMode(GL_MODELVIEW);
}

void GLText::drawFreeString(const char *string, float x, float y)
{
  	glColor3fv(fontColor);
   	glDisable(GL_DEPTH_TEST);
  	gl_font(fontType, fontSize );
  	gl_draw(string, x, y );
}

void GLText::drawFreeString(const char *string)
{
   	glColor3fv(fontColor);
   	glDisable(GL_DEPTH_TEST);
  	gl_font(fontType, fontSize );
  	gl_draw(string, 0.0f, 0.0f );
   	glEnable(GL_DEPTH_TEST);
}


