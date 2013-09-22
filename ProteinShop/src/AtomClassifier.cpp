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
AtomClassifier - Functions mapping atoms to a range of integers.

***********************************************************************/


#include <vector>
using namespace std;

#include "AtomClassifier.h"
using namespace MD;


/////////////////////////////////////////
/// AtomClassifier default singletons ///
/////////////////////////////////////////


class UnityClassifier : public AtomClassifier
{
public:

    ~UnityClassifier() {}
    uint classify (const Protein::ChainAtom &atom,
                   const ProteinState &proteinState)
    {
        return 0;
    }
    const char *className (uint classNum) { return "Unity"; }
    const char *name() { return "Unity"; }
    uint numClasses() { return 1; }
};
static UnityClassifier f_unitySingleton;


class PhobiPhilicClassifier : public AtomClassifier
{
public:

    ~PhobiPhilicClassifier() {}
    
    uint classify (const Protein::ChainAtom &atom,
                   const ProteinState &proteinState)
    {
        Protein::Residue *residue = atom.getResidue();
        Protein::Residue::ResidueType residueType = residue->getResType();
        if ( residueType == Protein::Residue::PHOBIC )
            return 0;
        else if ( residueType == Protein::Residue::PHILIC )
            return 2;
        else
            return 1;
    }

    const char *className (uint classNum)
    {
        static const char *classNames[3] = {
            "Hydrophobic",
            "Hydroneutral",
            "Hydrophilic"
        };
        if ( classNum < 3 ) return classNames[classNum];
        else return "<parameter range error>";
    }

    const char *name() { return "Phobic-Philic"; }
    uint numClasses() { return 3; }
};
static PhobiPhilicClassifier f_phobiPhilicSingleton;


class HBondClassifier : public AtomClassifier
{
public:

    ~HBondClassifier() {}

    uint classify (const MD::Protein::ChainAtom &atom,
                   const ProteinState &proteinState)
    {
        // first see if this atom is part of an amide group
        Protein::Residue *residue (atom.getResidue());
        Protein::Dipole dipole (residue->getAmide());
        Protein::ConstDipoleIterator dItr;
        if ( dipole.isValid() && (
                dipole.getMajorAtom() == &atom ||
                dipole.getMinorAtom() == &atom
             ) )
        {
            // got an amide, look for a carboxyl
            for ( dItr = proteinState.protein->carboxylsBegin();
                  dItr != proteinState.protein->carboxylsEnd();
                  ++dItr )
                if ( formHydrogenBond(dipole, *dItr) )
                    return 1;
            return 0;
        }
        // not an amide, see if it's in a carboxyl group
        dipole = residue->getCarboxyl();
        if ( dipole.isValid() && (
                dipole.getMajorAtom() == &atom ||
                dipole.getMinorAtom() == &atom
             ) )
        {
            // got a carboxyl, look for an amide
            for ( dItr = proteinState.protein->amidesBegin();
                  dItr != proteinState.protein->amidesEnd();
                  ++dItr )
                if ( formHydrogenBond(*dItr, dipole) )
                    return 2;
            return 0;
        }
        // not a dipole atom, no hydrogen bond
        return 0;
    }

    const char *className (uint classNum)
    {
        static const char *classNames[3] = {
            "Not Bonded",
            "Bonded Amide",
            "Bonded Carboxyl"
        };
        if ( classNum < 3 ) return classNames[classNum];
        else return "<parameter range error>";
    }

    const char *name() { return "Hydrogen Bond"; }
    uint numClasses() { return 3; }
};
static HBondClassifier f_hBondSingleton;


// copy and fill in this template to make a new classifier
/*
class Classifier : public AtomClassifier
{
public:

    ~Classifier() {}
    uint classify (const MD::Protein::ChainAtom &atom,
                   const ProteinState &proteinState)
    {
    }
    const char *className (uint classNum)
    {
    }
    const char *name()
    {
    }
    uint numClasses()
    {
    }
};
static Classifier f_Singleton;
*/


static const uint f_numDefaultClassifiers = 3;
static AtomClassifier *f_defaultClassifiers[f_numDefaultClassifiers] = {
    &f_unitySingleton,
    &f_phobiPhilicSingleton,
    &f_hBondSingleton
};
static vector<AtomClassifier*> f_userClassifiers;


//////////////////////////////
/// AtomClassifier methods ///
//////////////////////////////


AtomClassifier::~AtomClassifier()
{
    for ( uint i = 0; i < f_userClassifiers.size(); ++i )
    {
        if ( f_userClassifiers[i] == this )
        {
            f_userClassifiers.erase (f_userClassifiers.begin() + i);
            return;
        }
    }
}


uint AtomClassifier::add (AtomClassifier *classifier)
{
    if ( !classifier ) return UINT_MAX;
    for ( uint i = 0; i < f_numDefaultClassifiers; ++i )
        if ( classifier == f_defaultClassifiers[i] )
            return UINT_MAX;
    for ( uint i = 0; i < f_userClassifiers.size(); ++i )
        if ( classifier == f_userClassifiers[i] )
            return UINT_MAX;
    f_userClassifiers.push_back (classifier);
    return ( f_userClassifiers.size() - 1 );
}


AtomClassifier *AtomClassifier::get (uint number)
{
    if ( number < f_numDefaultClassifiers )
        return f_defaultClassifiers[number];
    number -= f_numDefaultClassifiers;
    if ( number < f_userClassifiers.size() )
        return f_userClassifiers[number];
    else
        return 0;
}


uint AtomClassifier::numClassifiers()
{
    return ( f_numDefaultClassifiers + f_userClassifiers.size() );
}


AtomClassifier *AtomClassifier::remove (uint number)
{
    if ( number < f_numDefaultClassifiers ) return 0;
    number -= f_numDefaultClassifiers;
    if ( number < f_userClassifiers.size() )
    {
        AtomClassifier *classifier = f_userClassifiers[number];
        f_userClassifiers.erase (f_userClassifiers.begin() + number);
        return classifier;
    }
    return 0;
}


/*
    Template fanatics may note that this code is identical to that found in
    ColorFunction.cpp -- these factories can be templatized for additional
    reuse, if another one is needed.
*/
