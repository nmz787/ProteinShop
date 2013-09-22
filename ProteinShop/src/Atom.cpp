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

#include <string.h>

#include "Atom.h"

namespace MD {

/*****************************
Static elements of class Atom:
*****************************/

const char Atom::elementNames[118][4]=
{"H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Uun","Uuu","Uub","UUt","UUq","Uup","Uuh","Uus","Uuo"};

const Scalar Atom::elementRadii[118]=
{
    0.53,0.31,
    1.67,1.12,0.87,0.67,0.56,0.48,0.42,0.38,
    1.90,1.45,1.18,1.11,0.98,0.88,0.79,0.71,
    2.43,1.94,1.84,1.76,1.71,1.66,1.61,1.56,1.52,1.49,1.45,1.42,1.36,1.25,1.14,1.03,0.94,0.88,
    2.65,2.19,2.12,2.06,1.98,1.90,1.83,1.78,1.73,1.69,1.65,1.61,1.56,1.45,1.33,1.23,1.15,1.08,
    2.98,2.53,0.00,0.00,2.47,2.06,2.05,2.38,2.31,2.33,2.25,2.28,2.26,2.26,2.22,2.22,2.17,2.08,2.00,1.93,1.88,1.85,1.80,1.77,1.74,1.71,1.56,1.54,1.43,1.35,1.27,1.20,
    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
};

const Scalar Atom::elementCovalentRadii[118]=
{
    0.37,0.32,
    1.34,0.90,0.82,0.77,0.75,0.73,0.71,0.69,
    1.54,1.30,1.18,1.11,1.06,1.02,0.99,0.97,
    1.96,1.74,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.19,1.16,1.14,1.10,
    2.11,1.92,1.62,1.48,1.37,1.45,1.56,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33,1.30,
    2.25,1.98,1.69,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.60,1.50,1.38,1.46,1.59,1.28,1.37,1.28,1.44,1.49,1.48,1.47,1.46,0.00,0.00,1.45,
    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
};

const Scalar Atom::elementVanDerWaalsRadii[118]=
{
    1.20,0.32,
    1.34,0.90,0.82,1.70,1.55,1.52,0.71,0.69,
    1.54,1.30,1.18,1.11,1.80,1.80,1.75,0.97,
    1.96,1.74,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.39,1.26,1.22,1.19,1.16,1.14,1.10,
    2.11,1.92,1.62,1.48,1.37,1.45,1.56,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33,1.30,
    2.25,1.98,1.69,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.60,1.50,1.38,1.46,1.59,1.28,1.37,1.28,1.44,1.49,1.48,1.47,1.46,0.00,0.00,1.45,
    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00
};

/*********************
Methods of class Atom:
*********************/

Atom::Element Atom::parseType(const char* elementName)
{
    /* Search through all element names: */
    int index;
    for(index=0;index<118;++index)
        if(strcmp(elementNames[index],elementName)==0)
            break;
    return Element(index);
}

void bond(Atom& atom1,Atom& atom2)
{
    /* Check if the bond already exists: */
    bool bondExists=false;
    for(std::vector<Atom*>::iterator bIt=atom1.bonds.begin();bIt!=atom1.bonds.end();++bIt)
        bondExists|=*bIt==&atom2;
    for(std::vector<Atom*>::iterator bIt=atom2.bonds.begin();bIt!=atom2.bonds.end();++bIt)
        bondExists|=*bIt==&atom1;
    
    if(!bondExists)
    {
        /* Create the bond: */
        atom1.bonds.push_back(&atom2);
        atom2.bonds.push_back(&atom1);
    }
}

}
