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

#include "Timer.h"

/**********************
Methods of class Timer:
**********************/

Timer::Timer(void)
    :elapsedSeconds(0),elapsedMicrons(0)
{
    struct timezone timeZone;
    gettimeofday(&lastMeasured,&timeZone);
}

void Timer::elapse(void)
{
    struct timeval newMeasured;
    struct timezone timeZone;
    gettimeofday(&newMeasured,&timeZone);
    
    elapsedMicrons=newMeasured.tv_usec-lastMeasured.tv_usec;
    int secondCarry=0;
    if(elapsedMicrons<0) // Correct for end-of-second wraparound
    {
        elapsedMicrons+=1000000;
        secondCarry=1;
    }
    elapsedSeconds=(newMeasured.tv_sec-secondCarry)-lastMeasured.tv_sec;
    if(elapsedSeconds<0) // Correct for end-of-day wraparound
        elapsedSeconds+=86400;
    lastMeasured=newMeasured;
}

double Timer::peekTime(void) const
{
    struct timeval newMeasured;
    struct timezone timeZone;
    gettimeofday(&newMeasured,&timeZone);
    
    int passedMicrons=newMeasured.tv_usec-lastMeasured.tv_usec;
    int secondCarry=0;
    if(passedMicrons<0) // Correct for end-of-second wraparound
    {
        passedMicrons+=1000000;
        secondCarry=1;
    }
    int passedSeconds=(newMeasured.tv_sec-secondCarry)-lastMeasured.tv_sec;
    if(passedSeconds<0) // Correct for end-of-day wraparound
        passedSeconds+=86400;
    
    return double(passedSeconds)+double(passedMicrons)/1000000.0;
}
