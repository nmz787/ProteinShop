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

#include <FL/fl_draw.H>

#include "Protein.h"
#include "ProteinRenderer.h"

#include "SequenceWidget.h"

/*******************************
Methods of class SequenceWidget:
*******************************/

void SequenceWidget::draw(void)
{
    if(protein!=0&&proteinRenderer!=0)
    {
        /* Draw residues as colored bars: */
        int ri1=protein->getResidueIndexRange().first;
        int ri2=protein->getResidueIndexRange().second;
        int lastX=x();
        for(int i=ri1;i<=ri2;++i)
        {
            int nextX=((i-ri1)*w())/(ri2-ri1)+x();
            const MD::Protein::Residue* rPtr=protein->pickResidue(i);
            MD::ProteinRenderer::Color rc=proteinRenderer->getResidueBackboneColor(rPtr);
            fl_color(uchar(rc[0]*255.0f*0.5f),uchar(rc[1]*255.0f*0.5f),uchar(rc[2]*255.0f*0.5f));
            fl_rectf(lastX,y(),nextX-lastX,h());
            lastX=nextX;
        }
    }
}
