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

#ifndef SERIALPORT_INCLUDED
#define SERIALPORT_INCLUDED

#include <string>
#include <utility>

class SerialPort
{
    /* Embedded classes: */
    public:
    enum ParitySettings
    {
        PARITY_NONE,PARITY_EVEN,PARITY_ODD
    };
    enum PortSettings
    {
        BLOCKING=0x0,NONBLOCKING=0x1
    };
    
    /* Elements: */
    private:
    int port; // I/O handle for the serial port
    bool initialized; // True if serial port was successfully initialized
    unsigned int totalBytesReceived,totalBytesSent; // Number of bytes sent/received
    unsigned int numReadSpins,numWriteSpins; // Number of unsuccessful read/write attempts
    
    /* Private methods: */
    void readBlocking(int numBytes,char* bytes); // Blocking read
    void writeBlocking(int numBytes,const char* bytes); // Blocking write
    
    /* Constructors and destructors: */
    public:
    SerialPort(std::string deviceName); // Opens the given device as a "raw" port for I/O
    ~SerialPort(void); // Closes the serial port
    
    /* Methods: */
    bool isInitialized(void) const // Checks for successful initialization
    {
        return initialized;
    };
    int getFd(void) const // Returns low-level file descriptor for serial port
    {
        return port;
    };
    void setPortSettings(int portSettingsMask); // Sets port file descriptor settings
    void setSerialSettings(int bitRate,int charLength,ParitySettings parity,int numStopbits,bool enableHandshake); // Sets serial port parameters
    void setRawMode(int minNumBytes,int timeout); // Switches port to "raw" mode and sets burst parameters
    void setCanonicalMode(void); // Switches port to canonical mode
    void setLineControl(bool respectModemLines,bool hangupOnClose); // Sets line control parameters
    unsigned int getTotalBytesReceived(void) const
    {
        return totalBytesReceived;
    };
    unsigned int getTotalBytesSent(void) const
    {
        return totalBytesSent;
    };
    unsigned int getNumReadSpins(void) const
    {
        return numReadSpins;
    };
    unsigned int getNumWriteSpins(void) const
    {
        return numWriteSpins;
    };
    void resetStatistics(void); // Resets the byte and spin counters
    bool waitForByte(double timeout); // Waits for a byte to appear on the serial port, returns true if byte can be read
    std::pair<char,bool> readByteNonBlocking(void); // Reads at most one byte from the port, second pair item is true if a byte was read
    char readByte(void) // Reads a single byte from the port
    {
        char result;
        readBlocking(1,&result);
        return result;
    };
    char* readBytes(int numBytes) // Reads a block of bytes; allocates a buffer using new[]
    {
        char* bytes=new char[numBytes];
        readBlocking(numBytes,bytes);
        return bytes;
    };
    char* readBytes(int numBytes,char* bytes) // Reads a block of bytes
    {
        readBlocking(numBytes,bytes);
        return bytes;
    };
    int readBytesRaw(int maxNumBytes,char* bytes); // Reads up to maximum number of available bytes; returns number of bytes read
    void writeByte(char byte) // Writes a single byte to the port
    {
        writeBlocking(1,&byte);
    };
    void writeBytes(int numBytes,const char* bytes) // Writes a block of bytes
    {
        writeBlocking(numBytes,bytes);
    };
    void writeString(std::string s) // Writes a standard string to the port
    {
        writeBlocking(s.length(),s.c_str());
    };
    void flush(void); // Flushes the output queue
    void drain(void); // Waits until all pending writes have completed
};

#endif
