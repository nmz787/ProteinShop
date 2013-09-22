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
Base Message Class.
***********************************************************************/

#ifndef _PSHOPMESSAGE_H_
#define _PSHOPMESSAGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* INFO */
enum PShopIDs{
    FIRSTID				= 100, 
    PROTEINID, 			/* Protein */
	RESIDUEID, 			/* Residue */
	AAID,				/* Amino Acid Sequence */
	SSID,				/* Second Structure	*/
	STRIDEID,			/* Stride Predicition */
	REDUCEID,			/* Reduce */
	PDBID,				/* PDB file I/O */
	PREDID,				/* Prediction file I/O */
	CONFIGID,			/* Configuration file I/O */
	ENGID,				/* Energy Function  */
	GEOMID,				/* Geometry */
	IKID,				/* Inverse Kinematic */
	NETID,				/* TCP client/server */
	SHUTDOWNID, 		/* Termination */
    INITID, 			/* Initialization */
    RUNID, 				/* Run time */
    /* Place new IDS before this line. */    
    LASTID
};

/* The base message from which all message classes must derive */
class PShopMessage {
public:
    /* Constructor */
    PShopMessage(void) {}; 
	virtual ~PShopMessage(void) {};
    int sendMessage( const char* kind, const int id, 
			    const char* stText, const char* exText, const char* argText);
	/* Fatal Message; program must terminate Leve 5 */
	int Fatal(const int id, const char* stText, const char* exText= NULL, 
			const char* argText= NULL) {
	sendMessage( "\n--Fatal", id, stText, exText, argText ); return 1; } 

	/* Error Message; error condition, can run without some features due to	error. Level 4*/
	int Error(const int id, const char* stText, const char* exText= NULL, 
			const char*	argText= NULL) {
	sendMessage("\n!!Error", id, stText, exText, argText ); return 1; } 
	
	/* Debug Message. Level 3 */
	int Debug(const int id, const char* stText, const char* exText= NULL, 
			const char* argText= NULL) {
	if (getenv("DEBUG")!=NULL)
	sendMessage( "\n~~Debug", id, stText, exText, argText); return 0; 
	}  

	/* Warning Message;  warn user of possible unexpected results. Level 2 */
	int Warn(const int id, const char* stText, const char* exText= NULL, 
			const char* argText= NULL) {
	sendMessage( "\n**Warning", id, stText, exText, argText );return 0;  }  

	/* Info Message; non-error information. Level 1 */
	int Info(const int id, const char* stText, const char* exText= NULL, 
			const char* argText= NULL) {
	sendMessage( "++Info", id, stText, exText, argText ); return 0; } 
    
	char buf[256];

private:
    char message[1024];
};


/* INLINE FUNCTION DEFINITIONS */
inline int PShopMessage::sendMessage( const char* kind, const int id, 
			    const char* stText, const char* exText, const char* argText)	{

    message[0] =  '\0';
	char cId[3];	// ids
    strcat(message, kind);
    sprintf(cId, " (%d)", id);
    strcat(message, cId);
    strcat(message, " : ");
    strcat(message, stText);
    if(exText)
	strcat(message, exText);
	if(argText)
	strcat(message, argText);
	strcat(message, "\n" );
	fprintf(stderr, "%s", message);
    fflush(stderr);
 	return 0;
} 

#endif // _PSHOPMESSAGE_H_
