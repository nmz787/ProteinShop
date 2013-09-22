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
ProteinFactory - Class to Create/Copy/Paste/Insert Protein .
***********************************************************************/

#include <cctype>

#include "ProteinFactory.h"
#include <FL/fl_ask.H>

namespace MD {

ProteinFactory::ProteinFactory(void)
{
	currentStructureType=Protein::SecondaryStructure::NONE;
}

ProteinFactory::~ProteinFactory(void)
{
}


Protein* ProteinFactory::copy(int beginIndex, int endIndex)
{
    ProteinState *state = curProtein();
    Protein *result = new Protein;
	Protein::ResidueCreator proteinCreator(result);
    Protein::Residue * residue;
	if ( !state || !state->protein->pickResidue(beginIndex)
				|| !state->protein->pickResidue(endIndex)) 
	{
		fl_message ("Please input valid residue indexes.");
		return 0;
	}
	for(int i = beginIndex; i < endIndex +1 ;i++)
	{
		residue = state->protein->pickResidue(i);
		factory(&proteinCreator, residue, i);
	}

    proteinCreator.finishProtein();
	return result;
}


void ProteinFactory::factory(ResidueCreator *proteinCreator, const Residue *residue,int index)
{
	if(currentStructureType != residue->getSecondaryStructure()->getStructureType())
	{
		currentStructureType =	residue->getSecondaryStructure()->getStructureType();
		switch (currentStructureType)
        	{
            	case Protein::SecondaryStructure::COIL:
         		  	proteinCreator->newSecondaryStructure(Protein::SecondaryStructure::COIL);
				break;

            	case Protein::SecondaryStructure::ALPHA_HELIX:
		  			proteinCreator->newSecondaryStructure(Protein::SecondaryStructure::ALPHA_HELIX);
				break;

            	case Protein::SecondaryStructure::BETA_STRAND:
		  			proteinCreator->newSecondaryStructure(Protein::SecondaryStructure::BETA_STRAND);
              	break;

	    		default:
	        	break;
        	}
		}
	proteinCreator->newResidue(residue->getPdbResidueName(),index);
	Geometry::Vector<double, 3> tmpPos;
	Position newPos;
  	Protein::ChainAtom* atomIt;
    for(atomIt=residue->beginAtom; atomIt!=residue->endAtom; atomIt=atomIt->succ)
    {
    	char elementName[2];
        char* atomNamePtr=atomIt->name;
        while(isspace(*atomNamePtr))
        	++atomNamePtr;
        elementName[0]=*atomNamePtr;
        ++atomNamePtr;
        elementName[1]='\0';

		proteinCreator->addAtom(elementName,atomIt->atomIndex,atomIt->getPosition(),atomNamePtr);
	}		

}

}
