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
MyRoller - Replacement for Fltk's roller widget with better callback
structure to allow undo buffers.
***********************************************************************/

#ifndef MYROLLER_INCLUDED
#define MYROLLER_INCLUDED

#include <FL/Fl_Roller.H>

class MyRoller:public Fl_Roller
{
    /* Methods inherited from Fl_Roller: */
    public:
    virtual int handle(int eventType);
    
    /* Constructors and destructors: */
    MyRoller(int x,int y,int w,int h,const char* label =0)
        :Fl_Roller(x,y,w,h,label)
    {
    };
};

#endif
