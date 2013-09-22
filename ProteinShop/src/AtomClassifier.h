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


class AtomClassifier;


#ifndef ATOM_CLASSIFIER_INCLUDED
#define ATOM_CLASSIFIER_INCLUDED


#include "Globals.h"
#include "Protein.h"


/** Interface for atom classifiers, which define arbitrary mappings from
    instances of MD::Protein::ChainAtom to a range of integers.  Some default
    implementations are provided, singleton instances of which can be obtained
    from the class factory.  In the future, the most useful implementations of
    this interface will be those that can be instantiated and modified by the
    end user (via a separate, as yet undefined, interface). */
class AtomClassifier
{
public:

    /** Remove this classifier from the factory's repertoire, if applicable. */
    virtual ~AtomClassifier();

    /** Classify an atom.
        @param atom
            The atom to classify.
        @param protein
            State data for the protein that @a atom belongs to.
        @return
            Integer in the range [0, numClasses()). */
    virtual uint classify (const MD::Protein::ChainAtom &atom,
                           const ProteinState &proteinState) = 0;

    /** Get the name of a class in this classifier.
        @param classNum
            Integer in the range [0, numClasses()).
        @return
            Name of the class suitable for user interfacing. */
    virtual const char *className (uint classNum) = 0;

    /** Get the name of this classifier.
        @return
            Name of this object suitable for user interfacing. */
    virtual const char *name() = 0;

    /** Get the range of this classifier.
        @return
            The number of classifications into which atoms are mapped. */
    virtual uint numClasses() = 0;

    /** Add a new classifier to the factory's repertoire.
        @param classifier
            Address of the object to add.  It will not be added if it is null or
            has already been added.
        @return
            The number that is assigned to this classifier, UINT_MAX if it is a
            duplicate. */
    static uint add (AtomClassifier *classifier);

    /** Factory method for atom classifiers.
        @param number
            Number of the classifier requested, in [0, numClassifiers()).
        @return
            Pointer to the requested classifier, or null if @a number is out of
            range. */
    static AtomClassifier *get (uint number);

    /** Get the number of classifiers available in the factory.
        @return
            The number of default classifiers, plus the number of classifiers
            that have been added, minus the number of classifiers that have been
            removed. */
    static uint numClassifiers();

    /** Remove a classifier from the factory's repertoire.
        @param number
            Number of the classifier requested, in [0, numClassifiers()).
        @return
            Pointer to the classifier that was removed, if @a number was in
            range; null otherwise.  The indices of all higher-numbered
            classifiers will decrease by 1 after removal. */
    static AtomClassifier *remove (uint number);
};


#endif
