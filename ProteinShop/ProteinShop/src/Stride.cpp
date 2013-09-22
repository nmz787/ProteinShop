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
Stride - Class to invoke the Stride 
***********************************************************************/

#define MAX_RESIDUE 4000

#include <stdexcept>
#include <stdio.h>      
#include <string.h>
#include <stdlib.h>
#include "Stride.h"
#include <typeinfo>
#include "Globals.h"

Stride::Stride(void)
	:count(0), selectChain(false)
{
     namebuffer= new char[MAX_RESIDUE];
	 namebuffer[0] ='\0';
     chainId[0] ='\0';
  	outFilename = ".stridePredFile";
}
 
Stride::~Stride(void)
{
	delete namebuffer;
}

int Stride::stridePrediction( const char *inFilename ) {

	#ifdef __APPLE__
  char *command   = "../plugin/Stride/stride.mac";
	#endif
	#ifdef __LINUX__
  char *command   = "../plugin/Stride/stride.linux";
	#endif
	#ifdef __SGI_IRIX__
  char *command   = "../plugin/Stride/stride.irix";
	#endif

  /* Check the availability of Stride binary */
  if (!command) {
	msg.Error(STRIDEID,"No Stride binary found!", command, "in ../plugin/Stride/bin" );
    return 0;
  }

  char *execuateStride = new char[strlen(command)+strlen(inFilename)+strlen(outFilename)+8];
  
  /* Generate command calls */
  sprintf(execuateStride,"%s %s -f%s", command, inFilename, outFilename);

  /* Validate the Stide prediction */
  if (system(execuateStride) <0 ){
  	msg.Fatal(STRIDEID,"Stride call failed");
    return 0;
  }
  
  /* Free memory */
  delete execuateStride;
  
  msg.Debug(STRIDEID,"Stride Second Structure prediction\n", namebuffer);

  return 1;
} 

 
//int Stride::parseStrideFile(const char *filename)
const char* Stride::parseStrideFile(const char *filename, const char* chainID )
{ 
	filename = ".stridePredFile";
  /* Set the specified chain for selection */
  if(chainID!=NULL) {
  	strncpy(chainId, chainID,1);
  	chainId[1] ='\0';
	msg.Debug(STRIDEID, chainId, " will be selected from the result of Stride.");
  }
  else
	msg.Debug(STRIDEID, "No chain selection for stride ", chainId);

	/* open the the input stride file and check for validity: */
	FILE *file = fopen(filename, "rt");
	if (!file) {
        msg.Error(STRIDEID,"Unable to open the prediction file	");
    	return NULL;
	}
  	
	char line[100], keyword[4], residuename[4], chain[1], origResIDS[5], ss[1], aa[50];
  	int  resID[1], currentID = 0;
	if(selectChain) 
  			msg.Debug(STRIDEID, chainId, " was specified for selection!");

	while (fgets(line, sizeof(line), file))
	{
         /* Parse the line just read: */
		sscanf(&line[ 0],"%3s",keyword);

		if (strcmp(keyword,"ASG")==0)
		{
			sscanf(line,"ASG  %3s %c %4s %4d    %c",
      			residuename, chain, origResIDS, resID, ss);
			
			/* Generate a sequence of Second Structure prediction from Stride */
			ss[1]='\0';
			if(!selectChain) 
			{
				strcat(namebuffer, ss);
				count++;
			}
			else
			{
				if(chainId != NULL) {
				/* Extract only the specified chain */
				if (strncmp(chain, chainId, 1)==0) 
				{
					strcat(namebuffer, ss);
					count++;
				}
				}
			}
		}
	}
    	
  	fclose(file);
	return namebuffer;
}  

void Stride::clean()
{
//  char *outfilename = ".stridePredFile";
  char *clean = new char[strlen("rm")+strlen(outFilename)+8];
  sprintf(clean,"\"%s\" %s", "rm", outFilename);
  /* Remove the output file from Stride */
  if (getenv("STRIDE")==NULL)
  system(clean);

  delete clean;

}

int Stride::writePredFile(const char *filename)
{ 
	/* open the the input stride file and check for validity: */
	FILE *file = fopen(filename, "w");
	if (!file) {
       // throw msg->Error(STRIDEID,"Stride: Unable to open input file", filename);
    	return 1;
	}
	fprintf(file,"Conf: ");
	for(int i = 0; i< count+1; i++) {
		if (i < count) {
		fprintf(file,"%c", '9');
		}
		else {
		fprintf(file,"\n");
		}
	}

  	fclose(file);
  	return 0;
}  

