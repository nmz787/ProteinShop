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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __APPLE__
#include <sys/types.h>
#endif

#include <sys/select.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdexcept>

#include "SpaceBall.h"

/**************************
Methods of class SpaceBall:
**************************/

bool SpaceBall::readLine(int lineBufferSize,char* lineBuffer,double timeout)
{
    int devicePortFd=devicePort.getFd();
    
    /* Initialize watched file descriptor sets: */
    fd_set readFds,writeFds,exceptFds;
    FD_ZERO(&readFds);
    FD_ZERO(&writeFds);
    FD_ZERO(&exceptFds);
    
    /* Compute "deadline" as current time + timeout: */
    struct timeval deadline;
    struct timezone timeZone;
    gettimeofday(&deadline,&timeZone);
    long timeoutSecs=long(floor(timeout));
    deadline.tv_sec+=timeoutSecs;
    deadline.tv_usec+=long(floor((timeout-double(timeoutSecs))*1000000.0));
    if(deadline.tv_usec>=1000000L)
    {
        ++deadline.tv_sec;
        deadline.tv_usec-=1000000L;
    }
    
    /* Return characters until carriage return is read or timeout expires: */
    int numRead=0;
    bool lineComplete=false;
    while(numRead<lineBufferSize-1)
    {
        /* Compute time interval until deadline: */
        struct timeval currentTime,timeoutTime;
        gettimeofday(&currentTime,&timeZone);
        timeoutTime.tv_sec=deadline.tv_sec-currentTime.tv_sec;
        timeoutTime.tv_usec=deadline.tv_usec-currentTime.tv_usec;
        if(timeoutTime.tv_usec<0L)
        {
            --timeoutTime.tv_sec;
            timeoutTime.tv_usec+=1000000L;
        }
        if(timeoutTime.tv_sec<0L) // Deadline expired?
            break;
        
        /* Wait until device ready for reading: */
        FD_SET(devicePortFd,&readFds);
        select(devicePortFd+1,&readFds,&writeFds,&exceptFds,&timeoutTime);
        if(FD_ISSET(devicePortFd,&readFds))
        {
            /* Read next available byte from the device port: */
            read(devicePortFd,lineBuffer+numRead,1);
            
            /* Check if we just read a CR or LF: */
            if(lineBuffer[numRead]=='\r'||lineBuffer[numRead]=='\n')
            {
                lineComplete=true;
                break;
            }
            else
                ++numRead;
        }
        else
            break;
    }
    
    /* Terminate read string and return: */
    lineBuffer[numRead]='\0';
    return lineComplete;
}

int SpaceBall::readPacket(int packetBufferSize,unsigned char* packetBuffer)
{
    /* Read characters until an end-of-line is encountered: */
    bool escape=false;
    int readBytes=0;
    while(readBytes<packetBufferSize-1)
    {
        /* Read next byte: */
        unsigned char byte=(unsigned char)(devicePort.readByte());

        /* Deal with escaped characters: */
        if(escape)
        {
            /* Process escaped character: */
            if(byte!='^') // Escaped circumflex stays
                byte&=0x1f;
            packetBuffer[readBytes]=byte;
            ++readBytes;
            escape=false;
        }
        else
        {
            /* Process normal character: */
            if(byte=='^') // Circumflex is escape character
                escape=true;
            else if(byte=='\r') // Carriage return denotes end of packet
                break;
            else
            {
                packetBuffer[readBytes]=byte;
                ++readBytes;
            }
        }
    }
    
    /* Terminate packet with ASCII NUL and return: */
    packetBuffer[readBytes]='\0'; 
    return readBytes;
}

void SpaceBall::deviceThread(void)
{
    /* Set thread's cancellation state: */
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,0);
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS,0);
    
    /* Receive lines from the serial port until interrupted: */
    while(true)
    {
        /* Read characters until an end-of-line is encountered: */
        unsigned char packet[256];
        int packetLength=readPacket(256,packet);
        
        /* Determine the packet type: */
        switch(packet[0])
        {
            case 'D':
            {
                /* Parse a data packet: */
                short int rawData[6];
                rawData[0]=(short int)(((unsigned int)packet[ 3]<<8)|(unsigned int)packet[ 4]);
                rawData[1]=(short int)(((unsigned int)packet[ 5]<<8)|(unsigned int)packet[ 6]);
                rawData[2]=(short int)(((unsigned int)packet[ 7]<<8)|(unsigned int)packet[ 8]);
                rawData[3]=(short int)(((unsigned int)packet[ 9]<<8)|(unsigned int)packet[10]);
                rawData[4]=(short int)(((unsigned int)packet[11]<<8)|(unsigned int)packet[12]);
                rawData[5]=(short int)(((unsigned int)packet[13]<<8)|(unsigned int)packet[14]);
                
                /* Convert axis states to floats: */
                float newAxisStates[6];
                for(int i=0;i<6;++i)
                {
                    bool neg=rawData[i]<0;
                    if(i%3==2)
                        neg=!neg;
                    newAxisStates[i]=float(abs(rawData[i]));
                    if(newAxisStates[i]>thresholds[i])
                        newAxisStates[i]-=thresholds[i];
                    else
                        newAxisStates[i]=0.0f;
                    newAxisStates[i]*=gains[i];
                    newAxisStates[i]=float(pow(newAxisStates[i],exponents[i]));
                    if(neg)
                        newAxisStates[i]=-newAxisStates[i];
                }
                
                /* Call event handler: */
                processAxisEvent(newAxisStates,axisStates);
                
                /* Update state arrays: */
                for(int i=0;i<6;++i)
                    axisStates[i]=newAxisStates[i];
                break;
            }
            
            case '.':
            {
                /* Parse a button event packet: */
                int buttonMask=0x0;
                buttonMask|=int(packet[2]&0x3f);
                buttonMask|=int(packet[2]&0x80)>>1;
                buttonMask|=int(packet[1]&0x1f)<<7;
                
                /* Convert button state mask to booleans: */
                bool newButtonStates[13];
                for(int i=0;i<13;++i)
                    newButtonStates[i]=buttonMask&(1<<i);
                
                /* Call event handler: */
                processButtonEvent(newButtonStates,buttonStates);
                
                /* Update state arrays: */
                for(int i=0;i<13;++i)
                    buttonStates[i]=newButtonStates[i];
                break;
            }
        }
    }
}

SpaceBall::SpaceBall(std::string devicePortName)
    :devicePort(devicePortName),
     numButtons(13),buttonStates(new bool[numButtons]),
     numAxes(6),
     thresholds(new float[numAxes]),gains(new float[numAxes]),exponents(new float[numAxes]),
     axisStates(new float[numAxes])
{
    /* Initialize gain and threshold arrays: */
    for(int i=0;i<numAxes;++i)
    {
        thresholds[i]=i>=3?10.0f:20.0f;
        gains[i]=i>=3?0.01f:0.001f;
        exponents[i]=2.0f;
    }
    
    /* Set device port parameters: */
    devicePort.setSerialSettings(9600,8,SerialPort::PARITY_NONE,2,false);
    devicePort.setRawMode(1,0);
    
    /* Wait for status message from device: */
    char lineBuffer[256];
    const int numResponses=4;
    char* responseTexts[numResponses]={"\021","@1 Spaceball alive and well","","@2 Firmware version"};
    int responseLengths[numResponses]={2,27,1,19};
    for(int i=0;i<numResponses;++i)
    {
        /* Try reading a line from the device port: */
        if(!readLine(256,lineBuffer,5.0))
            throw std::runtime_error("SpaceBall: Timeout while reading status message");
        
        /* Check if line is correct SpaceBall response: */
        if(strncmp(lineBuffer,responseTexts[i],responseLengths[i])!=0)
            throw std::runtime_error("SpaceBall: Incorrect response while reading status message");
    }
    
    /* Enable event reporting: */
    devicePort.writeString("M\r");
    
    /* Start device thread: */
    pthread_create(&deviceThreadId,0,deviceThreadWrapper,this);
}

SpaceBall::~SpaceBall(void)
{
    /* Stop device thread: */
    pthread_cancel(deviceThreadId);
    pthread_join(deviceThreadId,0);
    
    /* Delete state arrays: */
    delete[] buttonStates;
    delete[] thresholds;
    delete[] gains;
    delete[] exponents;
    delete[] axisStates;
}

void SpaceBall::processButtonEvent(const bool newButtonStates[],const bool oldButtonStates[])
{
}

void SpaceBall::processAxisEvent(const float newAxisStates[],const float oldAxisStates[])
{
}
