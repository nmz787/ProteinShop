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
SetDihedrals - Reads angles and residues from arrays and builds a PDB
protein structure.
***********************************************************************/

/* Change this to !=0 to write output file "AlphaBeta.pdb": */
#define WRITEPDBFILE 0

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "Protein.h"
#include "CreateProtein.h"
#define PI 3.14159265358979323

extern double residueAtomPos[25][25][3];
extern char residueAtomName[25][25][5];
extern double residued2[25];
extern double residued3[25];
extern double StandardPhi[25];
extern double StandardPsi[25];
extern double StandardAlpha[25];
extern double StandardBeta[25];
extern double StandardGamma[25];

extern int numbAtoms[25];

void identity4 (double a[4][4]) {
   int i, j;
   for (i = 0; i < 4; ++i) {
     for (j = 0; j < 4; ++j) a[i][j] = 0.;
     a[i][i] = 1.;
     }
   }

void translate4(double a[3], double b[4][4]) {
   int i;
   identity4(b);
   for (i = 0; i < 3; ++i) b[i][3] = a[i];
   }

void translate4d(double x, double y, double z, double b[4][4]) {
   int i;
   identity4(b);
   b[0][3] = x;
   b[1][3] = y;
   b[2][3] = z;
   }

void matmult4 (double a[4][4], double b[4][4], double c[4][4]) {

   int i, j, k;

   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 4; ++k)
            c[i][j] += a[i][k] * b[k][j];
         }
   }

void matmult_transp4 (double a[4][4], double b[4][4], double c[4][4]) {

   int i, j, k;

   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
         c[i][j] = 0.;
         for (k = 0; k < 4; ++k)
            c[i][j] += a[i][k] * b[j][k];
         }
   }

void matrix_vector4 (double a[4][4], double b[4], double c[4]) {

   int i, k;

   for (i = 0; i < 4; ++i) {
      c[i] = 0.;
      for (k = 0; k < 4; ++k)
         c[i] += a[i][k] * b[k];
      }
   }

void rotX4(double angle, double a[4][4]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
         a[i][j] = 0;
   a[1][1] = a[2][2] = cosa;
   a[1][2] = -sina;
   a[2][1] = sina;
   a[0][0] = 1.;
   a[3][3] = 1.;
   }

void rotY4(double angle, double a[4][4]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
         a[i][j] = 0;
   a[0][0] = a[2][2] = cosa;
   a[2][0] = -sina;
   a[0][2] = sina;
   a[1][1] = 1.;
   a[3][3] = 1.;
   }

void rotZ4(double angle, double a[4][4]) {

   int i, j;
   double sina, cosa;

   cosa = cos(angle);
   sina = sin(angle);
   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
         a[i][j] = 0;
   a[1][1] = a[0][0] = cosa;
   a[0][1] = -sina;
   a[1][0] = sina;
   a[2][2] = 1.;
   a[3][3] = 1.;
   }

namespace MD {

Protein* SetDihedrals (int numResidues, const int type[],
      const char pred[], const double phi[], const double psi[]) {
   int i, j, k, l, m, n;
   double v0[4], v1[4], v2[4], v3[4], v4[4], pos[4];
     double mainpos[4]={0.0,0.0,0.0,0.0};
   double a[4][4], b[4][4], c[4][4], d[4][4], e[4][4], f[4][4], g[4][4];
   #if WRITEPDBFILE
   FILE *fp;
   const char outfile[14] = "AlphaBeta.pdb";
   #endif
   Protein* result=new Protein;
   int currentResidueIndex=-1;
   Protein::ResidueCreator proteinCreator(result);
     Protein::SecondaryStructure::StructureType currentStructureType=Protein::SecondaryStructure::NONE;
   char chainId[2], pc;
   char elementName[2];
   int residueIndex;
   char* atomNamePtr;


   if(numResidues <= 0) return result;
   #if WRITEPDBFILE
   fp = fopen(outfile, "w");
   if (fp == 0) {
      printf("unable to open file %s\n", outfile);
      assert (fp != 0);
      }
   #endif
  n = 1;
  identity4(a);

//  Standard amino acids have CA at the origin and N at -residue_d2,
//  so translate to put N at the origin. Standard amino acids also
//  have C in the xy plane, but not necessarily carbonyl O.

    for (i = -1; i < numResidues+1; ++i)
    {
        if(i==-1||i==numResidues)
        {
            char* residuePdbName;
            if (i == -1)
            {
                j = 0;
                residuePdbName="ACE";
            }
            else
            {
                j = 14;
                residuePdbName="NME";
            }
            Protein::SecondaryStructure::StructureType newStructureType;
            newStructureType=Protein::SecondaryStructure::COIL;
            if(newStructureType!=currentStructureType)
            {
                proteinCreator.newSecondaryStructure(newStructureType);
                currentStructureType=newStructureType;
            }
            proteinCreator.newResidue(residuePdbName,i+1);
            for (k = 0; k <  numbAtoms[j]; ++k)
            {
                for (m = 0; m < 3; ++m)
                    pos[m] = residueAtomPos[j][k][m] + mainpos[m];
                pos[1] += 1;
                pos[0] -= 0.8;
                #if WRITEPDBFILE
                fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                        n, residueAtomName[j][k], residuePdbName,
                        i+1, pos[0], pos[1], pos[2], 0., 0.);
                #endif
                ++n;
                atomNamePtr=residueAtomName[j][k];
                while(isspace(*atomNamePtr))
                    ++atomNamePtr;
                elementName[0]=*atomNamePtr;
                ++atomNamePtr;
                elementName[1]='\0';
                proteinCreator.addAtom(elementName, n-1 ,Position(pos), atomNamePtr);
            }
        }
        else
        {
            j = type[i];
            pc = pred[i];
            Protein::SecondaryStructure::StructureType newStructureType;
            if (pc == 'C')
                newStructureType=Protein::SecondaryStructure::COIL;
            else if (pc == 'H')
                newStructureType=Protein::SecondaryStructure::ALPHA_HELIX;
            else if (pc == 'E')
                newStructureType=Protein::SecondaryStructure::BETA_STRAND;
            else
            {
                printf("Unknown secondary structure type.\n");
                exit(-1);
            }
            if(newStructureType!=currentStructureType)
            {
                proteinCreator.newSecondaryStructure(newStructureType);
                currentStructureType=newStructureType;
            }
            proteinCreator.newResidue(Protein::Residue::abbreviatedNames[j],i+1);

            translate4d(residued2[j], 0., 0., b);
            matmult4(a, b, c);
            rotX4(PI-StandardPhi[j], b);
            //   printf("StandardPhi[%d] = %f\n", j, StandardPhi[j]);
            matmult4(c, b, f);              // f is for the amide H

            //  First do NH

            if(j == 16)
                l = 1; 
            else
                l = 2;
            for (k = 0; k < l; ++k)
            {
                for (m = 0; m < 3; ++m)
                    v0[m] = residueAtomPos[j][k][m];
                v0[3] = 1.;
                if(k == 0) matrix_vector4(c, v0, mainpos);
                else  matrix_vector4(f, v0, mainpos);
                #if WRITEPDBFILE
                fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                        n, residueAtomName[j][k], Protein::Residue::abbreviatedNames[j],
                        i+1, mainpos[0], mainpos[1], mainpos[2], 0., 0.);
                #endif
                ++n;

                atomNamePtr=residueAtomName[j][k];
                while(isspace(*atomNamePtr))
                    ++atomNamePtr;
                elementName[0]=*atomNamePtr;
                ++atomNamePtr;
                elementName[1]='\0';
                proteinCreator.addAtom(elementName, n-1 ,Position(mainpos), atomNamePtr);
                //  printf("%d %4s %4s %4s\n",
                //            n-1, elementName, atomNamePtr, residueAtomName[j][k]);
            }
            rotX4(phi[i], b);
            matmult4(c, b, a);              // a continuing backbone
            rotZ4(StandardAlpha[j], b);
            matmult4(a, b, c);
            rotX4(psi[i], b);
            matmult4(c, b, d);              // d is to help continue backbone
            rotX4(StandardPsi[j], b);      // to bring O to plane of next N and CA
            matmult4(d, b, g);
            rotZ4(-StandardAlpha[j], b);
            matmult4(g, b, e);              // e is for carbonyl oxygen

            //  Now do rest of atoms

            for (k = l; k <  numbAtoms[j]; ++k)
            {
                for (m = 0; m < 3; ++m)
                    v0[m] = residueAtomPos[j][k][m];
                v0[3] = 1.;
                if (k == l+3)
                    matrix_vector4(e, v0, pos);
                else
                    matrix_vector4(a, v0, pos);
                #if WRITEPDBFILE
                fprintf(fp, "ATOM %6d %4s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                        n, residueAtomName[j][k], Protein::Residue::abbreviatedNames[j],
                        i+1, pos[0], pos[1], pos[2], 0., 0.);
                #endif
                ++n;
                atomNamePtr=residueAtomName[j][k];
                while(isspace(*atomNamePtr))
                    ++atomNamePtr;
                elementName[0]=*atomNamePtr;
                ++atomNamePtr;
                elementName[1]='\0';
                proteinCreator.addAtom(elementName, n-1 ,Position(pos), atomNamePtr);
            }
            // Now build up rest of main chain effect on matrix a

            translate4d(residued3[j], 0., 0., b);
            matmult4(d, b, c);
            rotZ4(StandardBeta[j], b);
            matmult4(c, b, d);
            translate4d(1.325, 0., 0., b);
            matmult4(d, b, c);
            rotX4(PI, b);
            matmult4(c, b, d);
            rotZ4(StandardGamma[j], b);
            matmult4(d, b, a);
            //   printf("Alpha %f  Beta %f  Gamma %f\n", 180*StandardAlpha[j]/PI,
            //      180*StandardBeta[j]/PI, 180*StandardGamma[j]/PI);
        }
    }

  #if WRITEPDBFILE
  fclose(fp);
  #endif
  proteinCreator.finishProtein();
  return result;
  }

}
