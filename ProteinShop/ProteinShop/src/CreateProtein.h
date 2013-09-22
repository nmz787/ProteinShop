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

#ifndef CREATEPROTEIN_INCLUDED
#define CREATEPROTEIN_INCLUDED

namespace MD {

/* Forward declarations: */
class Protein;

void ReadStandards(const char* standardsDirectory); // Reads standard configuration of amino acid residues
int  checkStandardsAtomNo(const char* standardsDirectory, const char* aminoName);
int  matchStandardsAtoms(const char* standardsDirectory, const char* aminoName, const char*atomName, int index);
void setAlphaHelixAngles(double phi,double psi); // Sets default dihedral angles for alpha helices in degrees
void setBetaStrandAngles(double phi,double psi); // Sets default dihedral angles for beta strands in degrees
void setCoilRegionAngles(double phi,double psi); // Sets default dihedral angles for coil regions in degrees
Protein* ReadPredictionFile(const char* predictionFilename, int setprotein); // Reads prediction file and returns protein structure
Protein* SetDihedrals(int numResidues,const int type[],const char pred[],const double phis[],const double psis[]); // Creates protein structure from given amino acid sequence and dihedral angles
Protein* loadProtein(const char* inputFileName); // Creates protein from input file; file type is determined by extension
}

#endif
