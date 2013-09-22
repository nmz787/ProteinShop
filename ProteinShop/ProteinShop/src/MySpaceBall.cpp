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

#include <Math/Math.h>

#include "MySpaceBall.h"

/****************************
Methods of class MySpaceBall:
****************************/

MySpaceBall::MySpaceBall(const ConfigurationFile::SectionIterator& configFileSection)
    :SpaceBall(configFileSection.retrieveString("./devicePort")),
     zTranslationMode(DOLLY),
     zTranslationModeButtonIndex(configFileSection.retrieveValue<int,ValueCoder<int> >("./zTranslationModeButtonIndex")),
     zTranslationModeButtonToggle(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./zTranslationModeButtonToggle")),
     dollyGain(configFileSection.retrieveValue<float,ValueCoder<float> >("./dollyGain")),
     maxOnlyButtonIndex(configFileSection.retrieveValue<int,ValueCoder<int> >("./maxOnlyButtonIndex")),
     maxOnlyButtonToggle(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./maxOnlyButtonToggle",false)),
     maxOnlyButtonState(false),
     modifiedAxisStates(new float[getNumAxes()]),
     dragButtonIndex(configFileSection.retrieveValue<int,ValueCoder<int> >("./dragButtonIndex")),
     dragButtonToggle(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./dragButtonToggle",false)),
     dragButtonState(false),
     axisEventPending(false),buttonEventPending(false)
{
    /* Adjust linear and angular gain factors: */
    float linearGain=configFileSection.retrieveValue<float,ValueCoder<float> >("./linearGain",1.0f);
    float angularGain=configFileSection.retrieveValue<float,ValueCoder<float> >("./angularGain",1.0f);
    for(int i=0;i<3;++i)
    {
        setGain(i,linearGain);
        setGain(3+i,angularGain);
    }
    
    /* Initialize modified axis states: */
    for(int i=0;i<getNumAxes();++i)
        modifiedAxisStates[i]=getAxisState(i);
}

MySpaceBall::~MySpaceBall(void)
{
    delete[] modifiedAxisStates;
}

void MySpaceBall::processButtonEvent(const bool newButtonStates[],const bool oldButtonStates[])
{
    /* Find the button that caused the event: */
    eventButtonIndex=-1;
    for(int i=0;i<getNumButtons();++i)
    {
        if(newButtonStates[i]!=oldButtonStates[i])
        {
            eventButtonIndex=i;
            break;
        }
    }
    
    /* Test for known button events: */
    if(eventButtonIndex==zTranslationModeButtonIndex)
    {
        if(zTranslationModeButtonToggle)
        {
            /* Toggle z translation mode on button release: */
            if(!newButtonStates[zTranslationModeButtonIndex])
            {
                if(zTranslationMode==PAN)
                    zTranslationMode=DOLLY;
                else
                    zTranslationMode=PAN;
            }
        }
        else
        {
            if(newButtonStates[zTranslationModeButtonIndex])
                zTranslationMode=PAN;
            else
                zTranslationMode=DOLLY;
        }
    }
    else if(eventButtonIndex==maxOnlyButtonIndex)
    {
        if(maxOnlyButtonToggle)
        {
            /* Toggle max axis button state on button release: */
            if(!newButtonStates[maxOnlyButtonIndex])
                maxOnlyButtonState=!maxOnlyButtonState;
        }
        else
            maxOnlyButtonState=newButtonStates[maxOnlyButtonIndex];
    }
    else if(eventButtonIndex==dragButtonIndex)
    {
        bool oldDragButtonState=dragButtonState;
        if(dragButtonToggle)
        {
            /* Toggle drag button state on button release: */
            if(!newButtonStates[dragButtonIndex])
                dragButtonState=!dragButtonState;
        }
        else
            dragButtonState=newButtonStates[dragButtonIndex];
        
        if(dragButtonState!=oldDragButtonState)
            buttonEventPending=true;
    }
    else if(eventButtonIndex>=0)
        buttonEventPending=true;
}

void MySpaceBall::processAxisEvent(const float newAxisStates[],const float oldAxisStates[])
{
    /* Save new axis states: */
    for(int i=0;i<getNumAxes();++i)
        modifiedAxisStates[i]=newAxisStates[i];
    
    if(maxOnlyButtonState)
    {
        /* Zero all axes but the maximum-value one: */
        int maxAxisIndex=-1;
        float maxAxisValue=0.0f;
        for(int i=0;i<getNumAxes();++i)
            if(Math::abs(modifiedAxisStates[i])>maxAxisValue)
            {
                maxAxisIndex=i;
                maxAxisValue=Math::abs(modifiedAxisStates[i]);
            }
        
        for(int i=0;i<getNumAxes();++i)
            if(i!=maxAxisIndex)
                modifiedAxisStates[i]=0.0f;
    }
    
    axisEventPending=true;
}
