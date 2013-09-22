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
ColorFunction - Map the unit interval to RGBA space.

***********************************************************************/


class ColorFunction;


#ifndef COLOR_FUNCTION_INCLUDED
#define COLOR_FUNCTION_INCLUDED


#include "GLTypes.h"


/** Function object for arbitrary mappings from the unit interval to RGBA
    space.  Some default implementations are provided, all but one of which set
    full opacity (A = 1).  Singleton instances of these can be obtained from the
    class factory. */
class ColorFunction
{
public:

    /** Mapping range type, component vector in RGBA space. */
    typedef GLColor<double,4> Color;

    /** Remove this function from the factory's repertoire, if applicable. */
    virtual ~ColorFunction();

    /** The color function.
        @param scalar
            Clamped to the interval [0,1] before mapping.
        @param result
            The color corresponding to @a scalar is stored here. */
    virtual void map (double scalar, Color &result) = 0;

    /** Provide a name for this function suitable for utilization by a user
        interface component.
        @return
            User-friendly name in a null-terminated string. */
    virtual const char *name() = 0;

    /** Add a new function to the factory's repertoire.
        @param function
            Address of the object to add.  It will not be added if it is null or
            has already been added.
        @return
            The number that is assigned to this function, UINT_MAX if it is a
            duplicate. */
    static unsigned int add (ColorFunction *function);

    /** Factory method for color functions.
        @param number
            Number of the function requested, in [0, numFunctions()).
        @return
            Pointer to the requested function, or null if @a number is out of
            range. */
    static ColorFunction *get (unsigned int number);

    /** Get the number of functions available in the factory.
        @return
            The number of default functions, plus the number of functions that
            have been added, minus the number of functions that have been
            removed. */
    static unsigned int numFunctions();

    /** Remove a function from the factory's repertoire.
        @param number
            Number of the function requested, in [0, numFunctions()).
        @return
            Pointer to the function that was removed, if @a number was in range;
            null otherwise.  The indices of all higher-numbered functions will
            decrease by 1 after removal. */
    static ColorFunction *remove (unsigned int number);
};


#endif
