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

#ifndef PROTEINAFCTORY_H
#define PROTEINFACTORY_H

#include "Protein.h"
#include "Globals.h"

namespace MD {

class  ProteinFactory : public Protein
{
	private:
	Protein::SecondaryStructure::StructureType currentStructureType;

	public:
	ProteinFactory(void);
	~ProteinFactory(void);
	Protein* copy(int beginIndex, int endIndex);
	void factory(ResidueCreator *proteinCreator, const Residue * residue, int index);
};

}
#endif
