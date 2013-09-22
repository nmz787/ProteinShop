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
Timer - A class to calculate CPU time usage across OSes.
***********************************************************************/

#ifndef TIMER_INCLUDED
#define TIMER_INCLUDED

#include <sys/time.h>

class Timer
{
    /* Elements: */
    private:
    struct timeval lastMeasured; // Time value at last measuring point
    int elapsedSeconds,elapsedMicrons; // Number of seconds and microseconds elapsed
    
    /* Constructors and destructors: */
    public:
    Timer(void); // Creates a timer and initializes with the current resource usage (starts timing)
    
    /* Methods: */
    void elapse(void); // Takes a snapshot of the current timer values
    int getSeconds(void) const // Returns the number of measured seconds
    {
        return elapsedSeconds;
    };
    int getMicrons(void) const // Returns the number of elapsed microseconds
    {
        return elapsedMicrons;
    };
    double getTime(void) const // Returns the amount of measured time, in seconds
    {
        return double(elapsedSeconds)+double(elapsedMicrons)/1000000.0;
    };
    double peekTime(void) const; // Returns the amount of time passed since the last time elapse() was called
};

#endif
