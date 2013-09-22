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
SerialPort - Class to ease programming serial ports under UNIX/Linux.
***********************************************************************/

#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>

#include "SerialPort.h"

/***************************
Methods of class SerialPort:
***************************/

void SerialPort::readBlocking(int numBytes,char* bytes)
{
    while(numBytes>0)
    {
        int bytesReceived=read(port,bytes,numBytes);
        numBytes-=bytesReceived;
        bytes+=bytesReceived;
        totalBytesReceived+=bytesReceived;
        if(numBytes>0)
            ++numReadSpins;
    }
}

void SerialPort::writeBlocking(int numBytes,const char* bytes)
{
    while(numBytes>0)
    {
        int bytesSent=write(port,bytes,numBytes);
        numBytes-=bytesSent;
        bytes+=bytesSent;
        totalBytesSent+=bytesSent;
        if(numBytes>0)
            ++numWriteSpins;
    }
}

SerialPort::SerialPort(std::string deviceName)
    :port(open(deviceName.c_str(),O_RDWR|O_NOCTTY)),initialized(port!=-1),
     totalBytesReceived(0),totalBytesSent(0),
     numReadSpins(0),numWriteSpins(0)
{
    if(initialized)
    {
        /* Configure as "raw" port: */
        struct termios term;
        tcgetattr(port,&term);
        term.c_iflag=IGNBRK|IGNPAR; // Don't generate signals or parity errors
        term.c_oflag=0; // No output processing
        term.c_cflag|=CREAD|CLOCAL; // Enable receiver; no modem line control
        term.c_lflag=0; // Don't generate signals or echos
        term.c_cc[VMIN]=1; // Block read() until at least a single byte is read
        term.c_cc[VTIME]=0; // No timeout on read()
        tcsetattr(port,TCSANOW,&term);
        
        /* Flush both queues: */
        tcflush(port,TCIFLUSH);
        tcflush(port,TCOFLUSH);
    }
}

SerialPort::~SerialPort(void)
{
    close(port);
}

void SerialPort::setPortSettings(int portSettingsMask)
{
    /* Retrieve current flags: */
    int fileFlags=fcntl(port,F_GETFL);
    
    /* Change flags according to given parameter: */
    if(portSettingsMask&NONBLOCKING)
        fileFlags|=FNDELAY|FNONBLOCK;
    else
        fileFlags&=~(FNDELAY|FNONBLOCK);
    
    /* Set new flags: */
    fcntl(port,F_SETFL,fileFlags);
}

void SerialPort::setSerialSettings(int bitRate,int charLength,SerialPort::ParitySettings parity,int numStopbits,bool enableHandshake)
{
    /* Initialize a termios structure: */
    struct termios term;
    tcgetattr(port,&term);
    
    /* Set rate of bits per second: */
    #ifdef __SGI_IRIX__
    term.c_ospeed=bitRate;
    term.c_ispeed=bitRate;
    #else
    /* Find the closest existing bitrate to the given one: */
    static speed_t bitRates[][2]={{0,B0},{50,B50},{75,B75},{110,B110},{134,B134},{150,B150},
                                  {200,B200},{300,B300},{600,B600},{1200,B1200},{1800,B1800},
                                  {2400,B2400},{4800,B4800},{9600,B9600},{19200,B19200},
                                  {38400,B38400},{57600,B57600},{115200,B115200},{230400,B230400}};
    int l=0;
    int r=19;
    while(r-l>1)
    {
        int m=(l+r)>>1;
        if(bitRate>=bitRates[m][0])
            l=m;
        else
            r=m;
    }
    cfsetospeed(&term,bitRates[l][1]);
    #endif
    
    /* Set character size: */
    term.c_cflag&=~CSIZE;
    switch(charLength)
    {
        case 5:
            term.c_cflag|=CS5;
            break;
        
        case 6:
            term.c_cflag|=CS6;
            break;
        
        case 7:
            term.c_cflag|=CS7;
            break;
        
        case 8:
            term.c_cflag|=CS8;
            break;
    }
    
    /* Set parity settings: */
    switch(parity)
    {
        case PARITY_ODD:
            term.c_cflag|=PARENB|PARODD;
            break;
        
        case PARITY_EVEN:
            term.c_cflag|=PARENB;
            break;
    }
    
    /* Set stop bit settings: */
    if(numStopbits==2)
        term.c_cflag|=CSTOPB;
    
    /* Set handshake settings: */
    if(enableHandshake)
    {
        #ifdef __SGI_IRIX__
        term.c_cflag|=CNEW_RTSCTS;
        #else
        term.c_cflag|=CRTSCTS;
        #endif
    }
        
    /* Set the port: */
    tcsetattr(port,TCSADRAIN,&term);
}

void SerialPort::setRawMode(int minNumBytes,int timeOut)
{
    /* Read the current port settings: */
    struct termios term;
    tcgetattr(port,&term);
    
    /* Disable canonical mode: */
    term.c_lflag&=~ICANON;
    
    /* Set the min/time parameters: */
    term.c_cc[VMIN]=cc_t(minNumBytes);
    term.c_cc[VTIME]=cc_t(timeOut);
    
    /* Set the port: */
    tcsetattr(port,TCSANOW,&term);
}

void SerialPort::setCanonicalMode(void)
{
    /* Read the current port settings: */
    struct termios term;
    tcgetattr(port,&term);
    
    /* Enable canonical mode: */
    term.c_lflag|=ICANON;
    
    /* Set the port: */
    tcsetattr(port,TCSANOW,&term);
}

void SerialPort::setLineControl(bool respectModemLines,bool hangupOnClose)
{
    /* Read the current port settings: */
    struct termios term;
    tcgetattr(port,&term);
    
    if(respectModemLines)
        term.c_cflag&=~CLOCAL;
    else
        term.c_cflag|=CLOCAL;
    if(hangupOnClose)
        term.c_cflag|=HUPCL;
    else
        term.c_cflag&=~HUPCL;
    
    /* Set the port: */
    tcsetattr(port,TCSANOW,&term);
}

void SerialPort::resetStatistics(void)
{
    totalBytesReceived=0;
    totalBytesSent=0;
    numReadSpins=0;
    numWriteSpins=0;
}

bool SerialPort::waitForByte(double timeout)
{
    /* Prepare parameters for select: */
    fd_set readFdSet;
    FD_ZERO(&readFdSet);
    FD_SET(port,&readFdSet);
    struct timeval tv;
    int seconds=int(floor(timeout));
    timeout-=double(seconds);
    int microseconds=int(floor(timeout*1000000.0+0.5));
    tv.tv_sec=seconds;
    tv.tv_usec=microseconds;
    
    /* Wait for an event on the port and return: */
    return select(port+1,&readFdSet,0,0,&tv)>0;
}

std::pair<char,bool> SerialPort::readByteNonBlocking(void)
{
    char readByte;
    int bytesReceived=read(port,&readByte,1);
    return std::pair<char,bool>(readByte,bytesReceived>0);
}

int SerialPort::readBytesRaw(int maxNumBytes,char* bytes)
{
    return read(port,bytes,maxNumBytes);
}

void SerialPort::flush(void)
{
    tcflush(port,TCOFLUSH);
}

void SerialPort::drain(void)
{
    tcdrain(port);
}
