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
ReadPredictionFile - Functions to read prediction files and create their
protein structures.
***********************************************************************/

#include <stdio.h>

#include "Protein.h"
#include "CreateProtein.h"
#include "Globals.h"
#include <stdexcept>

#define PI 3.14159265358979323

static int strand = 0, nstrand = 0;
static int coil = 0, ncoil = 0;
static int strandstart[100], strandend[100];
static int coilstart[100], coilend[100];

#define MAXRES 1400
static int type[MAXRES] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
char pred[MAXRES] = {'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'} ;
static int numResidues = 10;
static double phi[MAXRES];
static double psi[MAXRES];
//Nelson' change
int conf[MAXRES], is_core[MAXRES], has_core = 0;
//static int conf[MAXRES];
static double Cphi = -182.*PI/180;
static double Cpsi = -182.*PI/180;
static double Ephi = -120.*PI/180.;
static double Epsi = 140.*PI/180.;
static double Hphi = -60.*PI/180.;
static double Hpsi = -40.*PI/180.;
static double ProlinePhi = -60.*PI/180.;
extern char inputfilename[];

namespace MD {

Protein* ReadPredictionFile(const char* predictionFilename, int setprotein)
{
    int c,d,i,j;
	

    FILE* fp = fopen(predictionFilename, "r");
    strcpy (inputfilename, predictionFilename);
    if (!fp)
    {
		throw std::runtime_error(std::string("Unable to open ")+std::string(predictionFilename));
        return 0;
    }
	else
		msg.Info(PREDID, "Loading ", predictionFilename);
    
	numResidues = 0;
    while((c = fgetc(fp)) != ':')
        ;
    while((c = fgetc(fp)) == ' ')
        ;
    while(c >= '0' && c <= '9' && numResidues < MAXRES)
    {
        conf[numResidues] = c - '0';
        c = fgetc(fp);
        ++numResidues;
    }
    i = 0;
    while((c = fgetc(fp)) != ':')
        ;
    while((c = fgetc(fp)) == ' ')
        ;
    while((c == 'C' || c == 'E' || c == 'H') && i < MAXRES)
    {
        pred[i] = c;
        if(c == 'C')
        {
            phi[i] = Cphi;
            psi[i] = Cpsi;
            strand = 0;
            if(!coil)
            {
                coil = 1;
                ++ncoil;
                coilstart[ncoil] = i;
            }
            else
                coilend[ncoil] = i;
        }
        else if(c == 'E')
        {
            phi[i] = Ephi;
            psi[i] = Epsi;
            coil = 0;
            if(!strand)
            {
                strand = 1;
                ++nstrand;
                strandstart[nstrand] = i;
            }
            else
                strandend[nstrand] = i;
        }
        else if(c == 'H')
        {
            phi[i] = Hphi;
            psi[i] = Hpsi;
            coil = 0;
            strand = 0;
        }
        c = fgetc(fp);
        ++i;
    }
    if (i != numResidues)
    {
        printf("numbers of confidences %d and predictions %d do not agree\n", numResidues, i);
        return 0;
    }
    i = 0;
    while((c = fgetc(fp)) != ':')
        ;
    while((c = fgetc(fp)) == ' ')
        ;
    while(c > 64 && i < MAXRES)
    {
        for( j = 0; j < 25; ++j)
        {
            d = MD::Protein::Residue::singleLetterNames[j];
            if (c == d || c + ('a' - 'A') == d)
            {
                type[i] = j;
                //if(j == 18) phi[i] = ProlinePhi;
                break;
            }
        }
        if(j == 25)
        {
            printf("unknown single letter residue name %c\n", c);
            return 0;
        }
        c = fgetc(fp);
        ++i;
    }
    if (i != numResidues)
    {
        printf("numbers of confidences %d and single letter names %d do not agree\n", numResidues, i);
        return 0;
    }
    printf("setprotein in ReadPredictionFile is %d\n", setprotein);
    if(setprotein)
	return SetDihedrals( numResidues, type, pred, phi, psi);
    else
        return 0;
}

}
