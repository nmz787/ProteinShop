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
CreateProtein - Functions to create protein structures from amino acid
sequences and secondary structure prediction sequences.
***********************************************************************/

#include <string.h>
#include <stdexcept>

#include "Protein.h"
#include "ParsePdbFile.h"
#include "CreateProtein.h"

namespace MD {

Protein* loadProtein(const char* inputFileName)
{
    /* Find input file extension: */
    const char* extension;
    for(const char* cPtr=inputFileName;*cPtr!='\0';++cPtr)
        if(*cPtr=='.')
            extension=cPtr+1;
    
    /* Determine input file type based on extension: */
    if(strcasecmp(extension,"pdb")==0)
    {
        /* Load a PDB file: */
        return parsePdbFile(inputFileName);
    }
    else if(strcasecmp(extension,"pred")==0)
    {
        /* Create protein from prediction file: */
        return ReadPredictionFile(inputFileName, 1);
    }
    else
        throw std::runtime_error(std::string("loadProtein: Unknown file extension in input file ")+std::string(inputFileName));
}

}
