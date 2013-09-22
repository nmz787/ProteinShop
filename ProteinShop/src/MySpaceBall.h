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
MySpaceball - Class to connect a spaceball device to ProtoShop.
***********************************************************************/

#ifndef MYSPACEBALL_INCLUDED
#define MYSPACEBALL_INCLUDED

#include <ConfigurationFile.h>

#include "SpaceBall.h"

class MySpaceBall:public SpaceBall
{
    /* Embedded classes: */
    public:
    enum ZTranslationMode // Enumerated type for modes of Z translation
    {
        PAN,DOLLY
    };
    
    /* Elements: */
    private:
    ZTranslationMode zTranslationMode; // Current Z translation mode for navigation
    int zTranslationModeButtonIndex; // Index of button toggling Z translation mode
    bool zTranslationModeButtonToggle; // Flag if Z translation mode button acts as toggle
    float dollyGain; // Gain factor for dolly movement
    int maxOnlyButtonIndex; // Index of button used for max axis only reporting
    bool maxOnlyButtonToggle; // Flag if the max axis button acts as a toggle
    bool maxOnlyButtonState; // Current state of max only button
    float* modifiedAxisStates; // Array of modified axis values
    int dragButtonIndex; // Index of button used for dragging
    bool dragButtonToggle; // Flag if the drag button acts as a toggle
    bool dragButtonState; // Current state of drag button
    bool axisEventPending; // Flag if the SpaceBall has reported an axis event
    bool buttonEventPending; // Flag if the SpaceBall has reported a button event
    int eventButtonIndex; // Index of button that caused last event
    
    /* Constructors and destructors: */
    public:
    MySpaceBall(const ConfigurationFile::SectionIterator& configFileSection); // Constructs SpaceBall device by reading a configuration file section
    virtual ~MySpaceBall(void);
    
    /* Methods from SpaceBall: */
    virtual void processButtonEvent(const bool newButtonStates[],const bool oldButtonStates[]);
    virtual void processAxisEvent(const float newAxisStates[],const float oldAxisStates[]);
    
    /* New methods: */
    ZTranslationMode getZTranslationMode(void) const // Returns current Z translation mode
    {
        return zTranslationMode;
    };
    float getDollyGain(void) const // Returns dolly gain
    {
        return dollyGain;
    };
    int getDragButtonIndex(void) const // Returns index of drag button
    {
        return dragButtonIndex;
    };
    bool getDragState(void) const // Returns current dragging state
    {
        return dragButtonState;
    };
    bool isAxisEventPending(void) // Checks whether the SpaceBall has reported an axis event; resets event flag on-the-fly
    {
        bool result=axisEventPending;
        axisEventPending=false;
        return result;
    };
    bool isButtonEventPending(void) // Checks whether the SpaceBall has reported a button event
    {
        bool result=buttonEventPending;
        buttonEventPending=false;
        return result;
    };
    int getEventButtonIndex(void) const // Returns index of button causing button event
    {
        return eventButtonIndex;
    };
    float getModifiedAxisState(int axisIndex) const // Returns a modified axis state
    {
        return modifiedAxisStates[axisIndex];
    };
};

#endif
