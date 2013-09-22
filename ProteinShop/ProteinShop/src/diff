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
SpaceBall - Base class for device drivers for 6-DOF joysticks
(SpaceBall 4000FLX).
***********************************************************************/

#ifndef SPACEBALL_INCLUDED
#define SPACEBALL_INCLUDED

#include <pthread.h>
#include <string.h>
#include <SerialPort.h>

class SpaceBall
{
    /* Elements: */
    private:
    SerialPort devicePort; // Serial port the tracking device hardware is connected to
    pthread_t deviceThreadId; // Thread ID of device communication thread
    int numButtons; // Number of buttons on device (12)
    bool* buttonStates; // C-style array of button states
    int numAxes; // Number of axes on device (6)
    float* thresholds; // C-style array of threshold values for each axis
    float* gains; // C-style array of gain factors for each axis
    float* exponents; // C-style array of exponents for each axis
    float* axisStates; // C-style array of axis values
    
    /* Private methods: */
    bool readLine(int lineBufferSize,char* lineBuffer,double timeout); // Tries to read a line from serial port with timeout; returns true if complete line read
    int readPacket(int packetBufferSize,unsigned char* packetBuffer); // Reads a space ball status packet from the serial port; returns number of read characters
    static void* deviceThreadWrapper(void* threadArg) // Wrapper for device communication thread
    {
        static_cast<SpaceBall*>(threadArg)->deviceThread();
        return 0;
    };
    void deviceThread(void); // Device communication thread
    
    /* Constructors and destructors: */
    public:
    SpaceBall(std::string devicePortName); // Constructs device driver for specified device port; throws runtime error if device does not respond
    virtual ~SpaceBall(void); // Destroys device driver
    
    /* Methods: */
    int getNumButtons(void) const // Returns number of buttons on device
    {
        return numButtons;
    };
    bool getButtonState(int buttonIndex) const // Returns current state of a button
    {
        return buttonStates[buttonIndex];
    };
    int getNumAxes(void) const // Returns number of axes on device
    {
        return numAxes;
    };
    float getAxisState(int axisIndex) const // Returns current state of an axis
    {
        return axisStates[axisIndex];
    };
    void setThreshold(int axisIndex,float newThreshold)
    {
        thresholds[axisIndex]=newThreshold;
    };
    void setGain(int axisIndex,float newGain)
    {
        gains[axisIndex]=newGain;
    };
    void setExponent(int axisIndex,float newExponent)
    {
        exponents[axisIndex]=newExponent;
    };
    
    /* Methods called from device thread when a button/axis event is reported by the device: */
    virtual void processButtonEvent(const bool newButtonStates[],const bool oldButtonStates[]); // Processes a button event
    virtual void processAxisEvent(const float newAxisStates[],const float oldAxisStates[]); // Processes an axis event
};

#endif
