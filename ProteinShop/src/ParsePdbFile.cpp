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

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

#include "Protein.h"
#include "ParsePdbFile.h"
#include "Stride.h"
#include "Globals.h"
#include "CreateProtein.h"
extern char inputfilename[];
extern int conf[], is_core[], has_core;

namespace MD {

Protein* parsePdbFile (const char *filename, bool valid, const char* chain,const int modelId)
    {
	Protein* result=new Protein;
    Protein::ResidueCreator proteinCreator(result);
    Stride* prd = new Stride;
    bool reorder = false, padHydrogen = false, firstSecondStructure = false; 
    bool startModel=false, chainSelection=false;
    bool pdbFileZeroBased=false, pdbWithRemark=false, firstRun=false;
	const char* finalSSPred;
	int currentResidueIndex=-1;
    char line[1024];
    int aaIndex =0, atomIndexOffset = 0, residueIndex, core = 0;
	char residueName[4];
	double atomPosition[3], chainAtomPosition[3];
	char currentSecondStructure, predSecondStructure;
	bool firstAtom = true, validResidue  = false;

	try {
    	padHydrogen = configFile->retrieveValue<bool,ValueCoder<bool> >("/EnergyCalculator/addCaps", false );	
    	// don't cap ACE/NME if this pdb file is not valid for energy
		if (!valid) padHydrogen = false;
	}
    catch(TagNotFoundError error) {
    	padHydrogen = false;
    }

	// check a specific chain was selected or not
	if (!chain || strncmp(chain, "", 1) == 0 ) 
		chainSelection = false;
	else 
		chainSelection = true;

	 /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
	
	if(file==0)
		throw std::runtime_error(std::string("Unable to open ")+std::string(filename));
    
    strcpy(inputfilename, filename);

     while(fgets(line,sizeof(line),file)!=0)
    {
        /* Parse the line just read: */
        char keyword[7];
        int atomIndex;
        char atomName[5];
        char secondaryStructureType[20];
   	  	char chainId[1];
      	bool filter = false;
		int badResidue;
    	int modelIndex;
	   	chainId[0]='\0';
		sscanf(&line[ 0],"%6s",keyword);

		if(strcmp(keyword,"MODEL")==0)
        {
     		sscanf(&line[12],"%d",&modelIndex);
			if(modelIndex!=modelId)
				startModel = false;
			else 
				{
					startModel = true;
					sprintf(msg.buf, "The model ( %d ) was selected. ", modelId);
					msg.Info(PDBID, msg.buf);
				}
		}
		// skip other models
		if (!startModel && modelId > 0)
		continue;
		 
		if(strcmp(keyword,"REMARK")==0)
        {
            sscanf(&line[6],"%s",secondaryStructureType);

            /* Check if this is the beginning of a new secondary structure: */
            if(strcmp(secondaryStructureType,"ALPHA-HELIX")==0) {
                proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::ALPHA_HELIX);
               	pdbWithRemark = true;
			}
			else if(strcmp(secondaryStructureType,"BETA-STRAND")==0) {
                proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::BETA_STRAND);
               	pdbWithRemark = true;
			}
			else if(strcmp(secondaryStructureType,"COIL")==0) {
                proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::COIL);
        	   	pdbWithRemark = true;
			}
	        else if(strcmp(secondaryStructureType, "BEGINCORE")==0 ) {
		   printf("BEGINCORE detected\n");
		   core = 1;
		   has_core = 1;
		   }
	        else if(strcmp(secondaryStructureType, "ENDCORE")==0 ) {
		   printf("ENDCORE detected\n");
		   core = 0;
		   }
		}
		else if(strcmp(keyword,"ATOM")==0||
           strcmp(keyword,"atom")==0 )
        {
			sscanf(&line[ 6],"%d",&atomIndex);
            sscanf(&line[12],"%4s",atomName);
            sscanf(&line[17],"%3s",residueName);
			sscanf(&line[21],"%c",chainId);
      		sscanf(&line[22],"%d",&residueIndex);
            sscanf(&line[30],"%lg %lg %lg",&atomPosition[0],&atomPosition[1],&atomPosition[2]);

           /* Run a stride prediction when using a pdb file without REMARK */
			if (!firstRun && !pdbWithRemark ) {
				firstRun = true;
				/* Pad ACE residue at the beginning of the Amino Acid Sequence for the energy calucation*/	
	    		if (padHydrogen){
					if(strcmp(residueName, "ACE")!=0) {
						msg.Info(PDBID, "Padding ACE and NME residues for the energy computation ");
		  				/* The first second structure is created for ACE padding */
						proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::COIL);
						currentSecondStructure ='C';
						atomIndexOffset = proteinCreator.padHydrogenAtoms(residueIndex,
							result->getNumAtoms(), Position(atomPosition), 0);
					}
					/* adjust the atom index */
					if (atomIndex < atomIndexOffset) {
					atomIndex+=atomIndexOffset;
					sprintf(msg.buf, "atomIndexOffset %d !",atomIndexOffset);
					msg.Debug(PDBID, msg.buf);
					}
				}
				
				msg.Debug(PDBID, "Computing the Second Structure of ", filename);
				
				if (chainSelection)
				{
					msg.Info(PDBID, "The chain (", chain, ") was selected. ");
					finalSSPred = prd->parseStrideFile( filename, chain );
					//finalSSPred = prd->stridePrediction( filename, chain );
				}
				else
					finalSSPred = prd->parseStrideFile( filename, NULL ); 
					//finalSSPred = prd->stridePrediction( filename, NULL ); 
				  	//printf("parse %s\n", finalSSPred);
				
				if (finalSSPred==NULL)
				return 0;
				
			}

			/* Multiple chains */
			if (chainSelection) {
				/* Make selection when there is a match */
				if (strncmp(chainId, chain, 1 )==0) {
					filter = true;
					prd->enableChainSelection();
					// reserve the last atom position of selected chain for NME
					chainAtomPosition[0] = atomPosition[0];
					chainAtomPosition[1] = atomPosition[1];
					chainAtomPosition[2] = atomPosition[2];
					msg.Debug(PDBID, "This chain (",chain,") was selected. ");
				}
				else
					filter = false; 
			}
            
			#ifdef pdbFileZeroBased
			/* Check if the PDB file's residue indices are zero-based: */
            if(residueIndex==0)
                pdbFileZeroBased=true;
            if(!pdbFileZeroBased)
                --residueIndex; // Hack to get around 1-based residue indices for optimizer
			#endif
		
		if (filter || !chainSelection || pdbWithRemark) {
            /* Is this atom part of the current residue? */
            if(currentResidueIndex!=residueIndex)
            {
				is_core[residueIndex] = core;
				/* Create the Second Structure based on the prediction form	Stride */
				if (!pdbWithRemark) {
					
					if (strcmp(residueName, "ACE")!=0)
					{
					predSecondStructure = finalSSPred[aaIndex];
					//Convert the Second Structure type so that can be a whole for IK
					if (predSecondStructure == 'B' || predSecondStructure == 'T')
						predSecondStructure = 'C';
					if (predSecondStructure == 'G' )
						predSecondStructure = 'H';
				
				/* Check if this is the beginning of a new secondary structure: */
				if(currentSecondStructure != predSecondStructure &&
	       			predSecondStructure != '\0')
	      		  {
					currentSecondStructure = predSecondStructure;
					
					switch(currentSecondStructure)
					{
						//case 'G': //310Helix
						case 'H': //AlphaHelix
		  					proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::ALPHA_HELIX);
		  				break;

						case 'E': //Strand
		  					proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::BETA_STRAND);
		  				break;

						//case 'B': //Bridge	
						//case 'T': //Turn		
						case 'C': //Coil
		  						proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::COIL);
		  				break;
					
						default:
						break;
					}
				  }
					aaIndex++;
				  }
				  else {
					// The pdb file has ACE. 
					currentSecondStructure ='C';
					proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::COIL);
					}
						
				}
				/* Check the exsitance of backbone. ProteinShop can't handle negative residueIndex!? */
				if((strcmp(residueName, "ACE")!=0) && (strcmp(residueName,"NME")!=0)) {
					if(strcmp(atomName, "N")!=0 || residueIndex<0) {
						validResidue = false;
						sprintf(msg.buf, "%d %s is not a valid residue for ProteinShop!",residueIndex, residueName);
						msg.Warn(PDBID, msg.buf);
					}
					else
						validResidue = true;
				}
				else {
						//There are ACE and NME in pdb file
						validResidue = true; // assume ACE and NME are both valid				
					}
				/* Add the current residue to the protein only when this residue is valid: */
				if(validResidue) {
                	msg.Debug(PDBID, residueName);
					if(filter)
						proteinCreator.newResidue(residueName,residueIndex,chainId);
					else if (modelId > 0)
						proteinCreator.newResidue(residueName,residueIndex,0, modelId);
					else
						proteinCreator.newResidue(residueName,residueIndex);
				 }
				currentResidueIndex=residueIndex;
            }
           /* Parse the element name and structure position: */
            char elementName[2];
            char* atomNamePtr=atomName;
            while(isspace(*atomNamePtr))
                ++atomNamePtr;
            elementName[0]=*atomNamePtr;
            ++atomNamePtr;
            elementName[1]='\0';

            /* Add the atom to the residue only when this residue is valid: */
            if(validResidue)
			proteinCreator.addAtom(elementName,atomIndex,Position(atomPosition),atomNamePtr);
		}
	  }
    }
	
	/* Pad NME residue at the end of the Amino Acid Sequence for the energy calucation*/	
	/* Create a new coil second structure for NME if last one is not a coil */
	if (padHydrogen && strcmp(residueName, "NME")!=0) {
		if(!pdbWithRemark && currentSecondStructure!= 'C')
			proteinCreator.newSecondaryStructure(Protein::SecondaryStructure::COIL);
		if (chainSelection) 
			proteinCreator.padHydrogenAtoms(residueIndex,
				result->getNumAtoms(), Position(chainAtomPosition), 1);
		else
			proteinCreator.padHydrogenAtoms(residueIndex,
				result->getNumAtoms(), Position(atomPosition), 1);

		padHydrogen = false;
	}

    /* Finish constructing the protein: */
    proteinCreator.finishProtein();

    /* set conf[] to default values */
    for (int i = 0; i < residueIndex + 2; ++i)
       conf[i] = 5;

    /* Clean up and return the constructed protein: */
    fclose(file);
	delete prd;
    return result;
}


static bool writePdbFile(const Protein& protein,FILE *file,bool writeStructure)
{

    #if 1
    /* Write bogus energy value into file if no structure is requested: */
    if(!writeStructure)
        fprintf(file,"%f\n",0.0);
    #endif

    /* Start checking for secondary structure boundaries: */
    const Protein::SecondaryStructure* currentSecondaryStructure=0;

    /* Write the protein's atom in ascending index order: */
    for(Protein::ConstAtomIterator aIt=protein.atomsBegin();aIt!=protein.atomsEnd();++aIt)
    {
        const Protein::Residue* rPtr=aIt->getResidue();
        if(rPtr->getSecondaryStructure()!=currentSecondaryStructure)
        {
            if(writeStructure)
            {
                fprintf(file,"REMARK      ");
                switch(rPtr->getSecondaryStructure()->getStructureType())
                {
                    case Protein::SecondaryStructure::ALPHA_HELIX:
                        fprintf(file,"ALPHA-HELIX\n");
                        break;

                    case Protein::SecondaryStructure::BETA_STRAND:
                        fprintf(file,"BETA-STRAND\n");
                        break;

                    case Protein::SecondaryStructure::COIL:
                        fprintf(file,"COIL\n");
                        break;
		    
		    case Protein::SecondaryStructure::NONE:
                        break;
			
                }
            }
            currentSecondaryStructure=rPtr->getSecondaryStructure();
        }
        char atomName[5];
        strcpy(atomName,aIt->getElementName());
        strcat(atomName,aIt->getPlacement());
        Position p=aIt->getPosition();
        fprintf(file,"ATOM   %4d  %-4s%3s %1s %3d    ",aIt->getAtomIndex(),atomName,rPtr->getPdbResidueName(),"",rPtr->getPdbResidueIndex()); 
        fprintf(file,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",double(p[0]),double(p[1]),double(p[2]),1.0,0.0);
    }
    fclose(file);

    return true;
}


bool writePdbFile(const Protein& protein,int fd,bool writeStructure)
{
    FILE *file = fdopen (fd, "wt");
    if ( !file )
        return false;
    else
        return writePdbFile (protein, file, writeStructure);
}


bool writePdbFile(const Protein& protein,const char* filename,bool writeStructure)
{
    /* Open the input file and check for validity: */
    FILE* file = fopen (filename, "wt");
    if ( !file )
        return false;
    else
        return writePdbFile (protein, file, writeStructure);
}
bool checkPdbFile(const char* filename) {

	return checkHydrogen(filename);
}

bool checkHydrogen (const char *filename)
{ 
    bool firstLoop=true;
	int currentResidueIndex=-1;
    char line[1024];
	int aminoAtomNo = 0;
	char residueName[4];
	char currentResidueName[4];
    int match;
	 /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
	
	if(file==0)
		throw std::runtime_error(std::string("Unable to check ")+std::string(filename));
	else
		msg.Debug(PDBID, "Checking the ordering of atom in ", filename);
    
     while(fgets(line,sizeof(line),file)!=0)
    {
        /* Parse the line just read: */
        char keyword[7];
        int atomIndex;
        char atomName[5];
   	  	char chainId[1];
    	int residueIndex;
		double atomPosition[3];
	   	
		sscanf(&line[ 0],"%6s",keyword);
        
		if(strcmp(keyword,"ATOM")==0)
        {
			//sscanf(&line[ 6],"%d",&atomIndex);
            sscanf(&line[12],"%4s",atomName);
            sscanf(&line[17],"%3s",residueName);
      		sscanf(&line[22],"%d",&residueIndex);

			/* Is this atom part of the current residue? */
            if(firstLoop) {
                currentResidueIndex=residueIndex;
				firstLoop = false;
				strcpy(currentResidueName, residueName);
			}
           
		   if( currentResidueIndex!=residueIndex)
		   	break;
		   else
		   	aminoAtomNo++;
	  
	  		if(MD::matchStandardsAtoms("Standards/", currentResidueName,
					atomName, aminoAtomNo )){
				continue;
			}
			else {
				fclose(file);
				return false;
			}
	  }
    }
	
	/* Clean up : */
    fclose(file);
	return true;
}

bool checkBackBones (const char *filename)
{ 
    bool firstLoop=true;
	int  currentResidueIndex=-1;
    char line[1024];
	int  backbone = 0;
    
	 /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
	
	if(file==0)
		throw std::runtime_error(std::string("Unable to check ")+std::string(filename));
	else
		msg.Debug(PDBID, "Checking the existance of backbone ", filename);
    
     while(fgets(line,sizeof(line),file)!=0)
    {
        /* Parse the line just read: */
        char keyword[7];
        char atomName[5];
    	int residueIndex;
		sscanf(&line[ 0],"%6s",keyword);
        
		if(strcmp(keyword,"ATOM")==0)
        {
            sscanf(&line[12],"%4s",atomName);
      		sscanf(&line[22],"%d",&residueIndex);
 
			if(firstLoop) {
                currentResidueIndex=residueIndex;
				firstLoop = false;
			}

			/* Is this atom part of the current residue? */
		   if( currentResidueIndex!=residueIndex)
		   	break;
		   else
		   {
		   	if(strcmp(atomName,"N")==0)
			backbone++;
			
			if(strcmp(atomName,"CA")==0)
			backbone++;

			if(strcmp(atomName,"C")==0)
			backbone++;
			}
	  }
    }

	/* Clean up and return the result: */
    fclose(file);
    if (backbone < 3)
		return false;
	else
		return true;
}
bool checkBackBone (const char *atomName)
{ 
	if(strcmp(atomName,"N")!=0)
		return false;
	else
		return true;
}

char* parseChains (const char *filename)
{ 
    char line[1024];
	char* currentChain = new char [1];
    char* chainbuffer= new char[10];
  	char* chainId  = new char[1];
	chainbuffer[0] ='\0';
    currentChain[0] ='\0';
	int atomConect[5];    
    
	 /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
	
	if(file==0)
		throw std::runtime_error(std::string("Unable to check ")+std::string(filename));
	else
		msg.Debug(PDBID, "Parsing the chain list ", filename);
    
     while(fgets(line,sizeof(line),file)!=0)
    {
        /* Parse the line just read: */
        char keyword[7];
        char atomName[5];
    	int residueIndex;

		sscanf(&line[ 0],"%6s",keyword);
        
		if(strcmp(keyword,"ATOM")==0)
        {
			sscanf(&line[21],"%c",chainId);
			chainId[1] = '\0';
	//printf("%s\n", chainId);
			if (strncmp(chainId, currentChain, 1)!=0)
			{
			 	strncpy(currentChain, chainId, 1);
				currentChain[1] = '\0';
				strcat(chainbuffer, chainId);
				msg.Debug(PDBID, chainbuffer);
			 }
	  	}
    }
	/* Clean up and return the result: */
    fclose(file);
	delete [] chainId;
	delete [] currentChain;
	return chainbuffer;
}

int parseModels (const char *filename)
{ 
    char line[1024];
    int modelIndex = 0;
	 /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
	
	if(file==0)
		throw std::runtime_error(std::string("Unable to check ")+std::string(filename));
	else
		msg.Debug(PDBID, "Parsing the Model list ", filename);
    
     while(fgets(line,sizeof(line),file)!=0)
    {
        /* Parse the line just read: */
        char keyword[7];
    	//int modelIndex;
		sscanf(&line[ 0],"%6s",keyword);
        
		if(strcmp(keyword,"MODEL")==0)
        {
      		sscanf(&line[12],"%d",&modelIndex);
			//printf("model id %d\n", modelIndex);
	  	}
    }
	/* Clean up and return the result: */
    fclose(file);
	return modelIndex;
}

}   // end namespace MD

