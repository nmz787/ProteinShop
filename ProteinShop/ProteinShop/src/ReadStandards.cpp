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
ReadStandards - Functions to read residues from template files of
standard amino acids.
***********************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "Protein.h"
#include "CreateProtein.h"
#include "Globals.h"
#include <stdexcept>

#define PI 3.14159265358979323

// #include <Geometry/MatrixRotations.h>

// #include <Geometry/AffineTransformation.h>

double residueAtomPos[25][25][3];
char residueAtomName[25][25][5];
double residued2[25];
double residued3[25];
double StandardPhi[25];
double StandardPsi[25];
double StandardAlpha[25];
double StandardBeta[25];
double StandardGamma[25];
int numbAtoms[25];
double v1[3], v2[3], v3[3], v4[3], v5[3], n1[3], n2[3], n3[3], a1[3], a2[3];
double a[3][3], b[3][3], c[3][3], d[3][3], va[3];
double temp[25][3];
double ps, ph, d2, d3, t1, t2, t3, t4, t5;

double dot(double a[3], double b[3]) {
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
  }

void normalize( double a[3]) {
  double d;
  d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  a[0] /= d;
  a[1] /= d;
  a[2] /= d;
  }
  
void cross( double a[3], double b[3], double c[3]) {
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
   }

void matmult (double a[3][3], double b[3][3], double c[3][3]) {

   int i, j, k;

   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 3; ++k)
            c[i][j] += a[i][k] * b[k][j];
         }
   }

void matmult_transp (double a[3][3], double b[3][3], double c[3][3]) {

   int i, j, k;

   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 3; ++k)
            c[i][j] += a[i][k] * b[j][k];
         }
   }

void matrix_vector (double a[3][3], double b[3], double c[3]) {

   int i, k;

   for (i = 0; i < 3; ++i) {
      c[i] = 0.;
      for (k = 0; k < 3; ++k)
         c[i] += a[i][k] *b[k];
      }
   }

void rotX(double angle, double a[3][3]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         a[i][j] = 0;
   a[1][1] = a[2][2] = cosa;
   a[1][2] = -sina;
   a[2][1] = sina;
   a[0][0] = 1.;
   }

void rotY(double angle, double a[3][3]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         a[i][j] = 0;
   a[0][0] = a[2][2] = cosa;
   a[2][0] = -sina;
   a[0][2] = sina;
   a[1][1] = 1.;
   }

void rotZ(double angle, double a[3][3]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         a[i][j] = 0;
   a[1][1] = a[0][0] = cosa;
   a[0][1] = -sina;
   a[1][0] = sina;
   a[2][2] = 1.;
   }


namespace MD {

void ReadStandards(const char* standardsDirectory) 
    {
    int i, j, k, l, m, PenultimateRotamers = 1;
    char filename[256];
    char residueName[4], head[12];
    
    /* Build filename */
    for (i = 0; i < 23; ++i) {
    //if( i == 5 || i == 9 || i == 23) continue;
   if( i == 21) continue;
   strcpy(filename, standardsDirectory);
    char buildResidueName[4];
    if(PenultimateRotamers == 1)
       strcpy(buildResidueName,Protein::Residue::abbreviatedNames[i]);
    else
       strcpy(buildResidueName,Protein::Residue::abbreviatedLcNames[i]);
    strcat(filename, buildResidueName);
    if(PenultimateRotamers == 0) {
       if(i ==  5) strcat(filename, ".Thiol");
       }
    if(i == 9) strcat(filename, ".NDProtonated");
    if(PenultimateRotamers == 1)
       strcat(filename, ".pdb");
    else
       strcat(filename, ".new.Standard");
       
    /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
    if(file==0) {
		throw std::runtime_error(std::string("Unable to open ")+std::string(filename));
    }
		msg.Debug(PREDID, "Loading Standard ", filename);
    
	char line[1024];
        j = -1;
    while(fgets(line,sizeof(line),file)!=0)
    {
        /* Parse the line just read: */
        char keyword[7];
        int atomIndex;
        sscanf(&line[0],"%11s",head);
        if(strcmp("ATOM", head)) continue;
        ++j;
        sscanf(&line[12],"%4s",residueAtomName[i][j]);
        sscanf(&line[17],"%3s",residueName);
        sscanf(&line[30],"%lg",&residueAtomPos[i][j][0]);
        sscanf(&line[38],"%lg",&residueAtomPos[i][j][1]);
        sscanf(&line[46],"%lg",&residueAtomPos[i][j][2]);
        if (strcasecmp(residueName, buildResidueName)) {
            printf("Residue names do not agree: %s %s\n%d %s %s\n",
                residueName, buildResidueName,
                j, residueAtomName[i][j], filename);
        }
    }
    numbAtoms[i] = j+1;
    double d2 = 0., d3 = 0., e2, e3;
    StandardGamma[i] = 59.*PI/180.;
    if(i == 16) {   // proline
       l = 1;    // CA position 
       m = 11;   // use CG
       }
    else {
       l = 2;    // CA position
       m = 1;    // use HN
       }
    
    for (k = 0; k < 3; ++k) {

// v1 is H-N  vector, or CG-N vector for proline 
// v2 is N-CA vector; 
// v3 is CA-C vector; 
// v4 is C-O  vector; 

        v1[k] = residueAtomPos[i][m][k] - residueAtomPos[i][0][k];
        v2[k] = residueAtomPos[i][l][k] - residueAtomPos[i][0][k];
        v3[k] = residueAtomPos[i][l+2][k] - residueAtomPos[i][l][k];
        v4[k] = residueAtomPos[i][l+3][k] - residueAtomPos[i][l+2][k];
        d2 += v2[k]*v2[k];
        d3 += v3[k]*v3[k];
        v5[k] = residueAtomPos[i][0][k];
        va[k] = residueAtomPos[i][l][k];
    }
    d2 = sqrt(d2);
    d3 = sqrt(d3);
    residued2[i] = d2;
    residued3[i] = d3;
    if(0) printf("residueName %s d2 %f d3 %f\n", residueName, d2, d3);

// Compute dihedral angles

        cross(v1, v2, n1);
        cross(v2, v3, n2);
        cross(v3, v4, n3);
        normalize(n1);
        normalize(n2);
        normalize(n3);
        cross(n1, n2, a1);
        cross(n2, n3, a2);
        t1 = dot(n1, n2);
        t2 = dot(n2, n3);
        t3 = dot(v2, v3);
        t4 = dot(a1, v2);
        t5 = dot(a2, v3);
        ph = acos(t1);
        ps = acos(t2);
        StandardAlpha[i] = acos(t3/(d2*d3));
        
    StandardBeta[i] = 65.*PI/180.;
    
// For proteins, dihedral angles are clockwise when viewed along bond.

    if(t4 > 0.) ph = -ph;
    if(t5 > 0.) ps = -ps;

// ph and ps differ by 180 degrees from the backbone angles.
    
    StandardPhi[i] = ph + PI;
    StandardPsi[i] = ps + PI;

#if 0
    if(i == -1)
         cout << Protein::Residue::abbreviatedNames[i]  << " " << StandardPhi[i]
           *180/PI << " " << StandardPsi[i]*180/PI  << " t2  " << t2 <<
           "  t5  " << t5 << " Alpha " << StandardAlpha[i]*180/PI << endl;
#endif

    double angx, angy, angz;
/*
    angy = atan2(-v2[2], -v2[0]);
    angz = asin(-v2[1]/d2);
*/
    angy = atan2(v2[2], v2[0]);
    angz = asin(v2[1]/d2);

    rotY(angy, a);
    rotZ(-angz, b);
    matmult(b, a, c);
//  printf("%f  %f  v2 %f  %f  %f\n", angy, angz, v2[0], v2[1], v2[2]);

    for (j = 0; j < numbAtoms[i]; ++j) {
        for (k = 0; k < 3; ++k)
            residueAtomPos[i][j][k] -= va[k];
        matrix_vector(c, residueAtomPos[i][j], temp[j]);
        if (j < 0) printf("%d %4s %3s %f %f %f\n",
            j, residueAtomName[i][j], residueName,
            temp[j][0], temp[j][1], temp[j][2]);
}
    for (k = 0; k < 3; ++k) 
        v4[k] = temp[l+2][k] - temp[l][k];
    angx = atan2(v4[2], v4[1]);
    rotX(-angx, d);
    for (j = 0; j < numbAtoms[i]; ++j) {
        matrix_vector(d, temp[j], v3);
        for (k = 0; k < 3; ++k)
           residueAtomPos[i][j][k] = v3[k] - va[k];
        matrix_vector(d, temp[j], residueAtomPos[i][j]);
        if (i == -1) 
                  printf("ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
            j, residueAtomName[i][j], residueName,
            1+1, residueAtomPos[i][j][0], residueAtomPos[i][j][1],
            residueAtomPos[i][j][2], 0., 0.);
}

    /* Clean up and return the constructed protein: */
    fclose(file);
    fflush(stdout);
}
    }

int checkStandardsAtomNo(const char* standardsDirectory, const char* aminoName) 
    {
    int i, j, k, l;
    char filename[256];
    char residueName[4], head[12];
    char tmpName[4];
	Protein::Residue::AminoAcid type = Protein::Residue::parseType(aminoName);
	strcpy( tmpName, Protein::Residue::abbreviatedLcNames[type]);

    /* Build filename */
    for (i = 0; i < 23; ++i) {
    if( i == 21) continue;
    
	strcpy(filename, standardsDirectory);
    char buildResidueName[4];

    strcpy(buildResidueName,Protein::Residue::abbreviatedLcNames[i]);
	if (strcmp(tmpName, buildResidueName)!=0) {
		msg.Debug(PREDID, buildResidueName," != ", tmpName );
		continue;
	}
	else
	strcat(filename, buildResidueName);
    
	
	if(i == 9) strcat(filename, ".NDProtonated");
    if(i ==  5) strcat(filename, ".Thiol");
    strcat(filename, ".new.Standard");
       
    /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
    if(file==0) {
		throw std::runtime_error(std::string("Unable to open ")+std::string(filename));
	}
		msg.Debug(PREDID, "Loading Standard ", filename);
    
	char line[1024];
        j = 0;
    while(fgets(line,sizeof(line),file)!=0)
    {
        /* Parse the line just read: */
        char keyword[7];
        int atomIndex;
        sscanf(&line[0],"%11s",head);
        if(strcmp("ATOM", head)) continue;
        ++j;
        sscanf(&line[12],"%4s",residueAtomName[i][j]);
        sscanf(&line[17],"%3s",residueName);
        sscanf(&line[30],"%lg",&residueAtomPos[i][j][0]);
        sscanf(&line[38],"%lg",&residueAtomPos[i][j][1]);
        sscanf(&line[46],"%lg",&residueAtomPos[i][j][2]);
        if (strcasecmp(residueName, buildResidueName)) {
            printf("Residue names do not agree: %s %s\n%d %s %s\n",
                residueName, buildResidueName,
                j, residueAtomName[i][j], filename);
        }
    }
    
	return j;
   /* Clean up and return the constructed protein: */
    fclose(file);
    fflush(stdout);
 }
} //check

int matchStandardsAtoms(const char* standardsDirectory, const char* aminoName,
						const char*atomName, int index) 
    {
    int i, atomIndex;
    char filename[256], residueName[4], stdAtomName[5], buildResidueName[4];
    char  head[7], line[1024], tmpName[4];
	
	Protein::Residue::AminoAcid type = Protein::Residue::parseType(aminoName);
	strcpy( tmpName, Protein::Residue::abbreviatedLcNames[type]);

    /* Build filename */
    for (i = 0; i < 23; ++i) {
    	if( i == 21) continue;
    
	strcpy(filename, standardsDirectory);

    strcpy(buildResidueName,Protein::Residue::abbreviatedLcNames[i]);
	if (strcmp(tmpName, buildResidueName)!=0) {
		msg.Debug(PREDID, buildResidueName," != ", tmpName );
		continue;
	}
	else
		strcat(filename, buildResidueName);
 	
	if(i == 9) strcat(filename, ".NDProtonated");
    if(i ==  5) strcat(filename, ".Thiol");
    	strcat(filename, ".new.Standard");
       
    /* Open the input file and check for validity: */
    FILE* file=fopen(filename,"rt");
    if(file==0) {
		throw std::runtime_error(std::string("Unable to open ")+std::string(filename));
	}
		msg.Debug(PREDID, "Loading Standard ", filename);
    	while(fgets(line,sizeof(line),file)!=0)
    	{
        /* Parse the line just read: */
 		sscanf(&line[0],"%6s",head);
			
			if(strcmp("ATOM", head)==0) 
			{
 			sscanf(&line[6],"%d",&atomIndex);
			sscanf(&line[12],"%4s",stdAtomName);
        	sscanf(&line[17],"%3s",residueName);
        		
				if (atomIndex == index ) {
					if (strcmp(atomName, stdAtomName)==0) 
					{
						//printf("Pass %d std %s Target %s in aa %s\n", index, stdAtomName, atomName, residueName);
						return true;
					}
				}
				else
					continue;
			}
			 // atom
		} //while
    
   /* Clean up : */
    fclose(file);
    fflush(stdout);
	return false;
 } //for
} //match

} //MD
