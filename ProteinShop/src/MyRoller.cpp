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

#include "MyRoller.h"

#include "Globals.h"

/*************************
Methods of class MyRoller:
*************************/

int MyRoller::handle(int eventType)
{
    int result=0;
    ProteinState *state = curProtein();
    if ( !state ) return result;

    switch(eventType)
    {
        case FL_PUSH:
            /* Start interaction: */
            if(state->interactor->isValid())
                undoBuffer.startInteraction (state->interactor->getStructure(),
                                             state->interactor->getUpdateDirection());
            result=Fl_Roller::handle(eventType);
            break;

        case FL_RELEASE:
            result=Fl_Roller::handle(eventType);

            /* Finish interaction: */
            if(state->interactor->isValid())
            {
                undoBuffer.finishInteraction();
                updateGui();
            }
            break;

        default:
            result=Fl_Roller::handle(eventType);
    }

    return result;
}
