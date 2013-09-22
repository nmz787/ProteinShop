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


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>
using namespace std;

#include "ColorFunction.h"
#include "Globals.h"


////////////////////////////////////////
/// ColorFunction default singletons ///
////////////////////////////////////////


inline void makeIntensityPatternColor (
    uint first, uint second, uint third, double param,
    ColorFunction::Color &result
)
{
    assert ( param >= 0.0 && param <= 1.0 );
    assert ( first < 3 && second < 3 && third < 3 );

    result[first] = 1.0;
    if ( param <= 0.5 )
    {
        result[second] = param * 2.0;
        result[third] = 0.0;
    }
    else
    {
        result[second] = 1.0;
        result[third] = (param * 2.0) - 1.0;
    }
    result[3] = 1.0;
}


// add an invisible color function (black, alpha = 0)
// also update default color in EnergyRenderer (add one)
class InvisibleColorFunction : public ColorFunction
{
public:

    virtual ~InvisibleColorFunction() {}
    virtual const char *name() { return "Invisible"; }
    virtual void map (double scalar, Color &result)
    {
        result[0] = result[1] = result[2] = result[3] = 0.0;
    }
};
static InvisibleColorFunction f_invisibleSingleton;


class RainbowColorFunction : public ColorFunction
{
public:

    virtual ~RainbowColorFunction() {}
    virtual const char *name() { return "Rainbow Hues"; }
    virtual void map (double scalar, Color &result)
    {
        clamp (scalar, 0.0, 1.0);
        double radian = 2.0 * M_PI * scalar;
        if ( radian <= 2.0 * M_PI / 3.0 )
            result[0] = cos (0.75 * radian);
        else if ( radian >= 4.0 * M_PI / 3.0 )
            result[0] = cos (0.75 * (2.0 * M_PI - radian));
        else
            result[0] = 0.0;
        result[1] = sin (0.75 * radian);
        if ( result[1] < 0.0 )
            result[1] = 0.0;
        result[2] = sin (0.75 * (radian - 2.0 * M_PI / 3.0));
        if ( result[2] < 0.0 )
            result[2] = 0.0;
        result[3] = 1.0;
    }
};
static RainbowColorFunction f_rainbowSingleton;


class RedConstantColorFunction : public ColorFunction
{
public:

    virtual ~RedConstantColorFunction() {}
    virtual const char *name() { return "Red Constant"; }
    virtual void map (double scalar, Color &result)
    {
        result[0] = result[3] = 1.0;
        result[1] = result[2] = 0.0;
    }
};
static RedConstantColorFunction f_redConstantSingleton;


class GreenConstantColorFunction : public ColorFunction
{
public:

    virtual ~GreenConstantColorFunction() {}
    virtual const char *name() { return "Green Constant"; }
    virtual void map (double scalar, Color &result)
    {
        result[1] = result[3] = 1.0;
        result[0] = result[2] = 0.0;
    }
};
static GreenConstantColorFunction f_greenConstantSingleton;


class BlueConstantColorFunction : public ColorFunction
{
public:

    virtual ~BlueConstantColorFunction() {}
    virtual const char *name() { return "Blue Constant"; }
    virtual void map (double scalar, Color &result)
    {
        result[2] = result[3] = 1.0;
        result[0] = result[1] = 0.0;
    }
};
static BlueConstantColorFunction f_blueConstantSingleton;


class WhiteConstantColorFunction : public ColorFunction
{
public:

    virtual ~WhiteConstantColorFunction() {}
    virtual const char *name() { return "White Constant"; }
    virtual void map (double scalar, Color &result)
    {
        result[0] = result[1] = result[2] = result[3] = 1.0;
    }
};
static WhiteConstantColorFunction f_whiteConstantSingleton;


class OrangeIntensityColorFunction : public ColorFunction
{
public:

    virtual ~OrangeIntensityColorFunction() {}
    virtual const char *name() { return "Orange Intensity Pattern"; }
    virtual void map (double scalar, Color &result)
    {
        clamp (scalar, 0.0, 1.0);
        makeIntensityPatternColor (0, 1, 2, scalar, result);
    }
};
static OrangeIntensityColorFunction f_orangeIntensitySingleton;


class GreenIntensityColorFunction : public ColorFunction
{
public:

    virtual ~GreenIntensityColorFunction() {}
    virtual const char *name() { return "Green Intensity Pattern"; }
    virtual void map (double scalar, Color &result)
    {
        clamp (scalar, 0.0, 1.0);
        makeIntensityPatternColor (1, 2, 0, scalar, result);
    }
};
static GreenIntensityColorFunction f_greenIntensitySingleton;


class BlueIntensityColorFunction : public ColorFunction
{
public:

    virtual ~BlueIntensityColorFunction() {}
    virtual const char *name() { return "Blue Intensity Pattern"; }
    virtual void map (double scalar, Color &result)
    {
        clamp (scalar, 0.0, 1.0);
        makeIntensityPatternColor (2, 1, 0, scalar, result);
    }
};
static BlueIntensityColorFunction f_blueIntensitySingleton;


class MagentaPinkIntensityColorFunction : public ColorFunction
{
public:

    virtual ~MagentaPinkIntensityColorFunction() {}
    virtual const char *name()
    {
        return "Magenta-Pink Intensity Pattern";
    }
    virtual void map (double scalar, Color &result)
    {
        clamp (scalar, 0.0, 1.0);
        makeIntensityPatternColor (0, 2, 1, scalar, result);
    }
};
static MagentaPinkIntensityColorFunction f_magentaPinkIntensitySingleton;


class YellowGreenIntensityColorFunction : public ColorFunction
{
public:

    virtual ~YellowGreenIntensityColorFunction() {}
    virtual const char *name()
    {
        return "Yellow-Green Intensity Pattern";
    }
    virtual void map (double scalar, Color &result)
    {
        clamp (scalar, 0.0, 1.0);
        makeIntensityPatternColor (1, 0, 2, scalar, result);
    }
};
static YellowGreenIntensityColorFunction f_yellowGreenIntensitySingleton;


class MagentaBlueIntensityColorFunction : public ColorFunction
{
public:

    virtual ~MagentaBlueIntensityColorFunction() {}
    virtual const char *name()
    {
        return "Magenta-Blue Intensity Pattern";
    }
    virtual void map (double scalar, Color &result)
    {
        clamp (scalar, 0.0, 1.0);
        makeIntensityPatternColor (2, 0, 1, scalar, result);
    }
};
static MagentaBlueIntensityColorFunction f_magentaBlueIntensitySingleton;


// copy and fill in this template to make a new color function
/*
class ColorFunction : public ColorFunction
{
public:

    virtual ~ColorFunction() {}
    virtual const char *name()
    {
    }
    virtual void map (double scalar, Color &result)
    {
    }
};
static ColorFunction f_Singleton;
*/


static const uint f_numDefaultFunctions = 12;
static ColorFunction *f_defaultFunctions[f_numDefaultFunctions] = {
    &f_invisibleSingleton,
    &f_rainbowSingleton,
    &f_redConstantSingleton,
    &f_greenConstantSingleton,
    &f_blueConstantSingleton,
    &f_whiteConstantSingleton,
    &f_orangeIntensitySingleton,
    &f_greenIntensitySingleton,
    &f_blueIntensitySingleton,
    &f_magentaPinkIntensitySingleton,
    &f_magentaBlueIntensitySingleton,
    &f_yellowGreenIntensitySingleton
};
static vector<ColorFunction*> f_userFunctions;


/////////////////////////////
/// ColorFunction methods ///
/////////////////////////////


ColorFunction::~ColorFunction()
{
    for ( uint i = 0; i < f_userFunctions.size(); ++i )
    {
        if ( f_userFunctions[i] == this )
        {
            f_userFunctions.erase (f_userFunctions.begin() + i);
            return;
        }
    }
}


uint ColorFunction::add (ColorFunction *function)
{
    if ( !function ) return UINT_MAX;
    for ( uint i = 0; i < f_numDefaultFunctions; ++i )
        if ( function == f_defaultFunctions[i] )
            return UINT_MAX;
    for ( uint i = 0; i < f_userFunctions.size(); ++i )
        if ( function == f_userFunctions[i] )
            return UINT_MAX;
    f_userFunctions.push_back (function);
    return ( f_userFunctions.size() - 1 );
}


ColorFunction *ColorFunction::get (uint number)
{
    if ( number < f_numDefaultFunctions )
        return f_defaultFunctions[number];
    number -= f_numDefaultFunctions;
    if ( number < f_userFunctions.size() )
        return f_userFunctions[number];
    else
        return 0;
}


uint ColorFunction::numFunctions()
{
    return ( f_numDefaultFunctions + f_userFunctions.size() );
}


ColorFunction *ColorFunction::remove (uint number)
{
    if ( number < f_numDefaultFunctions ) return 0;
    number -= f_numDefaultFunctions;
    if ( number < f_userFunctions.size() )
    {
        ColorFunction *function = f_userFunctions[number];
        f_userFunctions.erase (f_userFunctions.begin() + number);
        return function;        
    }
    return 0;
}

