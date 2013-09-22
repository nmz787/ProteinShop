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

#ifndef _STRIDE_H
#define _STRIDE_H
 
class Stride {
 
public:
 
    Stride(void);
    ~Stride(void);
        
    int 	stridePrediction( const char *infilename );
    const char* 	parseStrideFile(const char *stridefile, const char* chainID); 
    int 	writePredFile(const char *stridefile); 
	void	enableChainSelection(void) { selectChain = true;}
	void 	clean(void);

private:
	char*	namebuffer;
	//char*	aaSeq;
	int 	count;
	bool   selectChain;
	char 	chainId[1];
  	char *outFilename;
};

#endif
