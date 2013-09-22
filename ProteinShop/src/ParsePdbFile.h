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
ParsePdbFile - Function to read a protein structure from a protein
database file.
***********************************************************************/

#ifndef PARSEPDBFILE_INCLUDED
#define PARSEPDBFILE_INCLUDED

namespace MD {

/* Forward declarations: */
class Protein;

Protein* parsePdbFile (const char* filename,bool valid=true, const char* chain=0, const int modelId =0);
char* parseChains (const char *filename);
int   parseModels (const char *filename);
bool writePdbFile(const Protein& protein,const char* filename,bool writeStructure =true);
bool writePdbFile(const Protein& protein,int fd,bool writeStructure =true);
bool checkPdbFile(const char* filename);
bool checkHydrogen (const char *filename);
bool checkBackBones (const char *filename);
bool checkBackBone (const char *atomName);
}

#endif
