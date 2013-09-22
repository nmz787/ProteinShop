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
SequenceWidget - Widget class to display and select residues and
secondary structures in sequence order.
***********************************************************************/

#ifndef SEQUENCEWIDGET_INCLUDED
#define SEQUENCEWIDGET_INCLUDED

#ifdef FLTK_NEW
#include <FL/Fl_Widget.H>
#else
#include <FL/Fl_Widget.h>
#endif

/* Forward declarations: */
namespace MD {
class Protein;
class ProteinRenderer;
}

class SequenceWidget:public Fl_Widget
{
    /* Elements: */
    private:
    MD::Protein* protein; // Protein whose sequence is displayed
    MD::ProteinRenderer* proteinRenderer; // Renderer for the same protein
    
    /* Methods inherited from Fl_Widget: */
    protected:
    virtual void draw(void);
    public:
    // virtual int handle(int eventType);
    
    /* Constructors and destructors: */
    public:
    SequenceWidget(int x,int y,int w,int h,const char* label =0)
        :Fl_Widget(x,y,w,h,label),
         protein(0),proteinRenderer(0)
    {
    };
    
    /* Methods: */
    void setProtein(MD::Protein* newProtein)
    {
        protein=newProtein;
    };
    void setProteinRenderer(MD::ProteinRenderer* newProteinRenderer)
    {
        proteinRenderer=newProteinRenderer;
    };
};

#endif
