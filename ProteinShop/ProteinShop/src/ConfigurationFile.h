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
ConfigurationFile - Class to handle permanent storage of configuration
data in human-readable text files.
***********************************************************************/

#ifndef CONFIGURATIONFILE_INCLUDED
#define CONFIGURATIONFILE_INCLUDED

/***********************************************************************
It seems that g++ version 3.2.0 has broken support for template method
overloading. Since the value retrieval methods with default value coder
parameters rely on that mechanism, we need to remove them for the broken
compiler version, alas.
***********************************************************************/

#ifdef __GNUC__
#define __CAN_OVERLOAD_TEMPLATES__ (__GNUC__ != 3 || __GNUC_MINOR__ != 2 || __GNUC_PATCHLEVEL__ != 0)
#else
#define __CAN_OVERLOAD_TEMPLATES__  1
#endif

#include <stdio.h>
#include <string.h>
#include <list>
#include <string>
#include <stdexcept>

/****************************************************************
Classes to encode/decode values into/from human-readable strings:
****************************************************************/

/* Exception class to report decoding errors: */
class DecodingError:public std::runtime_error
{
    /* Constructors and destructors: */
    public:
    DecodingError(const std::string& what_arg)
        :std::runtime_error(what_arg)
    {
    };
};

/* Generic value coder class: */
template <class ValueParam>
class ValueCoder
{
};

/* Specialized value coder classes for "standard" data types: */
template <>
class ValueCoder<std::string>
{
    /* Methods: */
    public:
    static std::string encode(const std::string& value)
    {
        return value;
    };
    static std::string decode(const std::string& text)
    {
        return text;
    };
};

template <>
class ValueCoder<bool>
{
    /* Methods: */
    public:
    static std::string encode(bool value)
    {
        if(value)
            return "true";
        else
            return "false";
    };
    static unsigned int decode(const std::string& text)
    {
        if(strcasecmp(text.c_str(),"true")==0)
            return true;
        else if (strcasecmp(text.c_str(),"false")==0)
            return false;
        else
            throw DecodingError(std::string("Unable to convert \"")+text+std::string("\" to bool"));
    };
};

template <>
class ValueCoder<unsigned int>
{
    /* Methods: */
    public:
    static std::string encode(unsigned int value)
    {
        char buffer[20];
        sprintf(buffer,"%u",value);
        return std::string(buffer);
    };
    static unsigned int decode(const std::string& text)
    {
        unsigned int result;
        if(sscanf(text.c_str(),"%u",&result)!=1)
            throw DecodingError(std::string("Unable to convert \"")+text+std::string("\" to unsigned int"));
        return result;
    };
};

template <>
class ValueCoder<int>
{
    /* Methods: */
    public:
    static std::string encode(int value)
    {
        char buffer[20];
        sprintf(buffer,"%d",value);
        return std::string(buffer);
    };
    static int decode(const std::string& text)
    {
        int result;
        if(sscanf(text.c_str(),"%d",&result)!=1)
            throw DecodingError(std::string("Unable to convert \"")+text+std::string("\" to int"));
        return result;
    };
};

template <>
class ValueCoder<float>
{
    /* Methods: */
    public:
    static std::string encode(float value)
    {
        char buffer[20];
        sprintf(buffer,"%g",value);
        return std::string(buffer);
    };
    static float decode(const std::string& text)
    {
        float result;
        if(sscanf(text.c_str(),"%f",&result)!=1)
            throw DecodingError(std::string("Unable to convert \"")+text+std::string("\" to float"));
        return result;
    };
};

template <>
class ValueCoder<double>
{
    /* Methods: */
    public:
    static std::string encode(double value)
    {
        char buffer[20];
        sprintf(buffer,"%g",value);
        return std::string(buffer);
    };
    static double decode(const std::string& text)
    {
        double result;
        if(sscanf(text.c_str(),"%lf",&result)!=1)
            throw DecodingError(std::string("Unable to convert \"")+text+std::string("\" to double"));
        return result;
    };
};

/* Exception class to report file IO errors: */
class FileOpenError:public std::runtime_error
{
    /* Constructors and destructors: */
    public:
    FileOpenError(const std::string& fileName)
        :std::runtime_error(std::string("Error opening file ")+fileName)
    {
    };
};

/* Exception class to report malformed configuration files: */
class MalformedConfigFileError:public std::runtime_error
{
    /* Constructors and destructors: */
    public:
    MalformedConfigFileError(const std::string& error,int lineNumber,const std::string& configFileName)
        :std::runtime_error(error+std::string("in line ")+ValueCoder<int>::encode(lineNumber)+std::string(" of file ")+configFileName)
    {
    };
};

/* Exception class to report missing tags: */
class TagNotFoundError:public std::runtime_error
{
    /* Constructors and destructors: */
    public:
    TagNotFoundError(const std::string& tag)
        :std::runtime_error(std::string("Tag \"")+tag+std::string("\" not found"))
    {
    };
};

class ConfigurationFile
{
    /* Embedded classes: */
    private:
    class TagValue // Class representing tag/value pairs
    {
        /* Elements: */
        public:
        std::string tag;
        std::string value; // Value encoded as std::string
        
        /* Constructors and destructors: */
        TagValue(const std::string& sTag,const std::string& sValue) // Creates a std::string value
            :tag(sTag),value(sValue)
        {
        };
    };
    
    class Section // Class representing a section of configuration data
    {
        /* Elements: */
        public:
        Section* parent; // Pointer to parent section (null if root section)
        std::string name; // Section name
        Section* sibling; // Pointer to next section under common parent
        Section* firstSubsection; // Pointer to first subsection
        Section* lastSubsection; // Pointer to last subsection
        std::list<TagValue> values; // List of values in this section
        
        /* Constructors and destructors: */
        Section(Section* sParent,const std::string& sName); // Creates an empty section
        ~Section(void);
        
        /* Methods: */
        Section* addSubsection(const std::string& subsectionName); // Adds a subsection to a section
        void save(FILE* file,int sectionLevel) const; // Writes all subsections and tag/value pairs to a file
    };
    
    public:
    class SectionIterator // Class allowing shortcut access to file sections without using globally current section
    {
        friend class ConfigurationFile;
        
        /* Elements: */
        private:
        ConfigurationFile* file; // File containing represented section
        Section* baseSection; // Base section for resolving relative tag names
        
        /* Constructors and destructors: */
        public:
        SectionIterator(void) // Constructs illegal section iterator
            :file(0),baseSection(0)
        {
        };
        private:
        SectionIterator(ConfigurationFile* sFile,Section* sBaseSection) // Constructs iterator pointing to given section in given file
            :file(sFile),baseSection(sBaseSection)
        {
        };
        
        /* Methods: */
        public:
        
        /* Section iterator methods: */
        friend bool operator==(const SectionIterator& sIt1,const SectionIterator& sIt2) // Equality operator
        {
            return sIt1.file==sIt2.file&&sIt1.baseSection==sIt2.baseSection;
        };
        friend bool operator!=(const SectionIterator& sIt1,const SectionIterator& sIt2) // Inequality operator
        {
            return sIt1.file!=sIt2.file||sIt1.baseSection!=sIt2.baseSection;
        };
        SectionIterator& operator++(void) // Moves to next subsection under same parent section (pre-increment)
        {
            baseSection=baseSection->sibling;
            return *this;
        };
        SectionIterator operator++(int) // Ditto (post-increment)
        {
            SectionIterator result=*this;
            baseSection=baseSection->sibling;
            return result;
        };
        SectionIterator beginSubsections(void) const // Returns iterator to first subsection
        {
            return SectionIterator(file,baseSection->firstSubsection);
        };
        SectionIterator endSubsections(void) const // Returns iterator one past last subsection
        {
            return SectionIterator(file,0);
        };
        
        /* Section navigation methods: */
        std::string getSection(void) const; // Returns represented section as absolute path
        void setSection(const char* newRelativePath); // Sets the represented section to the given section; does its best to create sections that do not already exist
        
        /* Retrieval methods that throw exceptions if a requested tag does not exist: */
        template <class ValueParam,class ValueCoderParam>
        ValueParam retrieveValue(const char* tag) const // Retrieves a value; throws exception if tag does not exist
        {
            return ValueCoderParam::decode(file->retrieveTagValue(baseSection,tag));
        };
        #if __CAN_OVERLOAD_TEMPLATES__
        template <class ValueParam>
        ValueParam retrieveValue(const char* tag) const // Ditto
        {
            return ValueCoder<ValueParam>::decode(file->retrieveTagValue(baseSection,tag));
        };
        #endif
        const std::string& retrieveString(const char* tag) const // Retrieves string value; throws exception if tag does not exist
        {
            return file->retrieveTagValue(baseSection,tag);
        };
        
        /* Retrieval methods that return default value (and add the tag) if a requested tag does not exist: */
        template <class ValueParam,class ValueCoderParam>
        ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) const // Retrieves a value; returns default if tag does not exist
        {
            return ValueCoderParam::decode(file->retrieveTagValue(baseSection,tag,ValueCoderParam::encode(defaultValue)));
        };
        #if __CAN_OVERLOAD_TEMPLATES__
        template <class ValueParam>
        ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) const // Ditto
        {
            return ValueCoder<ValueParam>::decode(file->retrieveTagValue(baseSection,tag,ValueCoder<ValueParam>::encode(defaultValue)));
        };
        #endif
        template <class ValueParam,class ValueCoderParam>
        ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) // Retrieves a value; returns default (and adds tag) if tag does not exist
        {
            return ValueCoderParam::decode(file->retrieveTagValue(baseSection,tag,ValueCoderParam::encode(defaultValue)));
        };
        #if __CAN_OVERLOAD_TEMPLATES__
        template <class ValueParam>
        ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) // Ditto
        {
            return ValueCoder<ValueParam>::decode(file->retrieveTagValue(baseSection,tag,ValueCoder<ValueParam>::encode(defaultValue)));
        };
        #endif
        std::string retrieveString(const char* tag,const std::string& defaultValue) const // Retrieves a string value; returns default if tag does not exist
        {
            return file->retrieveTagValue(tag,defaultValue);
        };
        std::string& retrieveString(const char* tag,const std::string& defaultValue) // Retrieves a string value; returns default (and adds tag) if tag does not exist
        {
            return file->retrieveTagValue(baseSection,tag,defaultValue);
        };
        
        /* Methods to store values: */
        template <class ValueParam,class ValueCoderParam>
        void storeValue(const char* tag,const ValueParam& newValue) // Stores a value; adds tag if tag does not exist
        {
            std::string newString=ValueCoderParam::encode(newValue);
            file->retrieveTagValue(baseSection,tag,newString)=newString; // A wee bit of double effort; if tag does not exist, it will be set twice. Who cares?
        };
        #if __CAN_OVERLOAD_TEMPLATES__
        template <class ValueParam>
        void storeValue(const char* tag,const ValueParam& newValue) // Ditto
        {
            std::string newString=ValueCoder<ValueParam>::encode(newValue);
            file->retrieveTagValue(baseSection,tag,newString)=newString; // A wee bit of double effort; if tag does not exist, it will be set twice. Who cares?
        };
        #endif
        void storeString(const char* tag,const std::string& newValue) // Stores a string value; adds tag if tag does not exist
        {
            file->retrieveTagValue(baseSection,tag,newValue)=newValue; // A wee bit of double effort; if tag does not exist, it will be set twice. Who cares?
        };
    };
    
    friend class SectionIterator;
    
    /* Elements: */
    private:
    std::string fileName; // File name of configuration file
    Section rootSection; // Root section of configuration file
    Section* currentSection; // Pointer to currently active section
    bool isEdited; // Flag indicating that the current in-memory version of the configuration file differs from permanent version
    
    /* Private methods: */
    Section* getSection(Section* baseSection,const char* newRelativePath); // Returns section reached from base section following given path
    const std::string& retrieveTagValue(const Section* baseSection,const char* relativeTagPath) const; // Retrieves value of relative tag path; throws exception if tag does not exist
    const std::string& retrieveTagValue(const char* relativeTagPath) const // Ditto starting from current section
    {
        return retrieveTagValue(currentSection,relativeTagPath);
    };
    std::string retrieveTagValue(const Section* baseSection,const char* relativeTagPath,const std::string& defaultValue) const; // Retrieves value of relative tag path; returns default value if tag does not already exist
    std::string retrieveTagValue(const char* relativeTagPath,const std::string& defaultValue) const // Ditto starting from current section
    {
        return retrieveTagValue(currentSection,relativeTagPath,defaultValue);
    };
    std::string& retrieveTagValue(Section* baseSection,const char* relativeTagPath,const std::string& defaultValue); // Retrieves value of relative tag path; does its best to create tag if it does not already exist
    std::string& retrieveTagValue(const char* relativeTagPath,const std::string& defaultValue) // Ditto starting from current section
    {
        return retrieveTagValue(currentSection,relativeTagPath,defaultValue);
    };
    
    /* Constructors and destructors: */
    public:
    ConfigurationFile(const char* sFileName); // Opens an existing configuration file
    ~ConfigurationFile(void);
    
    /* Methods: */
    void load(void); // Reloads contents of configuration file
    void save(void); // Saves the current in-memory state of the configuration file
    
    /* Current section management methods: */
    std::string getCurrentSection(void) const; // Returns current section
    void setCurrentSection(const char* newCurrentSection); // Sets the current section to the given section; does its best to create sections that do not already exist
    void list(void) const; // Lists all subsections and tags in current section
    
    /* Section iterator management methods: */
    SectionIterator getRootSectionIterator(void) // Returns root section as section iterator
    {
        return SectionIterator(this,&rootSection);
    };
    SectionIterator getCurrentSectionIterator(void) // Returns current section as section iterator
    {
        return SectionIterator(this,currentSection);
    };
    SectionIterator getSectionIterator(const char* relativePath) // Returns section iterator for given relative path
    {
        return SectionIterator(this,getSection(currentSection,relativePath));
    };
    void setCurrentSection(const SectionIterator& newCurrentSection) // Sets current section to section iterator
    {
        if(newCurrentSection.file==this) // Only set section if it belongs to same file!
            currentSection=newCurrentSection.baseSection;
        else
            throw std::runtime_error("Attempt to set current section from different configuration file");
    };
    
    /* Retrieval methods that throw exceptions if a requested tag does not exist: */
    template <class ValueParam,class ValueCoderParam>
    ValueParam retrieveValue(const char* tag) const // Retrieves a value; throws exception if tag does not exist
    {
        return ValueCoderParam::decode(retrieveTagValue(tag));
    };
    #if __CAN_OVERLOAD_TEMPLATES__
    template <class ValueParam>
    ValueParam retrieveValue(const char* tag) const // Ditto
    {
        return ValueCoder<ValueParam>::decode(retrieveTagValue(tag));
    };
    #endif
    const std::string& retrieveString(const char* tag) const // Retrieves string value; throws exception if tag does not exist
    {
        return retrieveTagValue(tag);
    };
    
    /* Retrieval methods that return default value (and add the tag) if a requested tag does not exist: */
    template <class ValueParam,class ValueCoderParam>
    ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) const // Retrieves a value; returns default if tag does not exist
    {
        return ValueCoderParam::decode(retrieveTagValue(tag,ValueCoderParam::encode(defaultValue)));
    };
    #if __CAN_OVERLOAD_TEMPLATES__
    template <class ValueParam>
    ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) const // Ditto
    {
        return ValueCoder<ValueParam>::decode(retrieveTagValue(tag,ValueCoder<ValueParam>::encode(defaultValue)));
    };
    #endif
    template <class ValueParam,class ValueCoderParam>
    ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) // Retrieves a value; returns default (and adds tag) if tag does not exist
    {
        return ValueCoderParam::decode(retrieveTagValue(tag,ValueCoderParam::encode(defaultValue)));
    };
    #if __CAN_OVERLOAD_TEMPLATES__
    template <class ValueParam>
    ValueParam retrieveValue(const char* tag,const ValueParam& defaultValue) // Ditto
    {
        return ValueCoder<ValueParam>::decode(retrieveTagValue(tag,ValueCoder<ValueParam>::encode(defaultValue)));
    };
    #endif
    std::string retrieveString(const char* tag,const std::string& defaultValue) const // Retrieves a string value; returns default if tag does not exist
    {
        return retrieveTagValue(tag,defaultValue);
    };
    std::string& retrieveString(const char* tag,const std::string& defaultValue) // Retrieves a string value; returns default (and adds tag) if tag does not exist
    {
        return retrieveTagValue(tag,defaultValue);
    };
    
    /* Methods to store values: */
    template <class ValueParam,class ValueCoderParam>
    void storeValue(const char* tag,const ValueParam& newValue) // Stores a value; adds tag if tag does not exist
    {
        std::string newString=ValueCoderParam::encode(newValue);
        retrieveTagValue(tag,newString)=newString; // A wee bit of double effort; if tag does not exist, it will be set twice. Who cares?
    };
    #if __CAN_OVERLOAD_TEMPLATES__
    template <class ValueParam>
    void storeValue(const char* tag,const ValueParam& newValue) // Ditto
    {
        std::string newString=ValueCoder<ValueParam>::encode(newValue);
        retrieveTagValue(tag,newString)=newString; // A wee bit of double effort; if tag does not exist, it will be set twice. Who cares?
    };
    #endif
    void storeString(const char* tag,const std::string& newValue) // Stores a string value; adds tag if tag does not exist
    {
        retrieveTagValue(tag,newValue)=newValue; // A wee bit of double effort; if tag does not exist, it will be set twice. Who cares?
    };
};

/************************************
Shortcuts for retrieve/store methods:
************************************/

template <class ValueParam>
inline ValueParam retrieveValue(const ConfigurationFile& configFile,const char* tag)
{
    return configFile.template retrieveValue<ValueParam,ValueCoder<ValueParam> >(tag);
}

inline const std::string& retrieveString(const ConfigurationFile& configFile,const char* tag)
{
    return configFile.retrieveString(tag);
}

template <class ValueParam>
inline ValueParam retrieveValue(const ConfigurationFile& configFile,const char* tag,const ValueParam& defaultValue)
{
    return configFile.template retrieveValue<ValueParam,ValueCoder<ValueParam> >(tag,defaultValue);
}

inline std::string retrieveString(const ConfigurationFile& configFile,const char* tag,const std::string& defaultValue)
{
    return configFile.retrieveString(tag,defaultValue);
}

template <class ValueParam>
inline ValueParam retrieveValue(ConfigurationFile& configFile,const char* tag,const ValueParam& defaultValue)
{
    return configFile.template retrieveValue<ValueParam,ValueCoder<ValueParam> >(tag,defaultValue);
}

inline std::string& retrieveString(ConfigurationFile& configFile,const char* tag,const std::string& defaultValue)
{
    return configFile.retrieveString(tag,defaultValue);
}

template <class ValueParam>
inline void storeValue(ConfigurationFile& configFile,const char* tag,const ValueParam& newValue)
{
    configFile.template storeValue<ValueParam>(tag,newValue);
}

inline void storeString(ConfigurationFile& configFile,const char* tag,const std::string& newValue)
{
    configFile.storeString(tag,newValue);
}

#endif
