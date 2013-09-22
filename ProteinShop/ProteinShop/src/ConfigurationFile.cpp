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

#include <ctype.h>

#include "ConfigurationFile.h"

/*******************************************
Methods of class ConfigurationFile::Section:
*******************************************/

ConfigurationFile::Section::Section(ConfigurationFile::Section* sParent,const std::string& sName)
    :parent(sParent),name(sName),sibling(0),firstSubsection(0),lastSubsection(0)
{
}

ConfigurationFile::Section::~Section(void)
{
    /* Delete all subsections: */
    while(firstSubsection!=0)
    {
        Section* next=firstSubsection->sibling;
        delete firstSubsection;
        firstSubsection=next;
    }
}

ConfigurationFile::Section* ConfigurationFile::Section::addSubsection(const std::string& subsectionName)
{
    /* Check if the subsection already exists: */
    Section* sPtr;
    for(sPtr=firstSubsection;sPtr!=0;sPtr=sPtr->sibling)
        if(sPtr->name==subsectionName)
            break;
    
    if(sPtr==0)
    {
        /* Add new subsection: */
        Section* newSubsection=new Section(this,subsectionName);
        if(lastSubsection!=0)
            lastSubsection=lastSubsection->sibling=newSubsection;
        else
            firstSubsection=lastSubsection=newSubsection;
        
        return newSubsection;
    }
    else
        return sPtr;
}

void ConfigurationFile::Section::save(FILE* file,int sectionLevel) const
{
    /* Generate indentation: */
    char prefix[80];
    for(int i=0;i<sectionLevel;++i)
        prefix[i]='\t';
    prefix[sectionLevel]='\0';
    
    /* Write subsections followed by tag/value pairs: */
    bool didWriteSomething=false;
    
    /* Write all subsections: */
    for(const Section* ssPtr=firstSubsection;ssPtr!=0;ssPtr=ssPtr->sibling)
    {
        /* Write separator line: */
        if(didWriteSomething)
            fprintf(file,"%s\n",prefix);
        
        /* Write section header: */
        fprintf(file,"%ssection %s\n",prefix,ssPtr->name.c_str());
        
        /* Write section contents: */
        ssPtr->save(file,sectionLevel+1);
        
        /* Write section footer: */
        fprintf(file,"%sendsection\n",prefix);
        
        didWriteSomething=true;
    }
    
    /* Write tag/value pairs: */
    for(std::list<TagValue>::const_iterator tvIt=values.begin();tvIt!=values.end();++tvIt)
    {
        /* Write separator line: */
        if(didWriteSomething)
            fprintf(file,"%s\n",prefix);
        
        /* Write tag/value pair: */
        fprintf(file,"%s%s %s\n",prefix,tvIt->tag.c_str(),tvIt->value.c_str());
    }
}

/***************************************************
Methods of class ConfigurationFile::SectionIterator:
***************************************************/

std::string ConfigurationFile::SectionIterator::getSection(void) const
{
    /* Assemble section prefix string on the way up from the current section: */
    std::string result;
    for(Section* sPtr=baseSection;sPtr!=0;sPtr=sPtr->parent)
        result=sPtr->name+std::string("/")+result;
    return result;
}

void ConfigurationFile::SectionIterator::setSection(const char* newRelativePath)
{
    baseSection=file->getSection(baseSection,newRelativePath);
}

/**********************************
Methods of class ConfigurationFile:
**********************************/

ConfigurationFile::Section* ConfigurationFile::getSection(ConfigurationFile::Section* baseSection,const char* newRelativePath)
{
    /* Determine the section to start the tag search from: */
    Section* sPtr;
    const char* pathSuffixPtr=newRelativePath;
    if(pathSuffixPtr[0]=='/')
    {
        sPtr=&rootSection;
        ++pathSuffixPtr;
    }
    else
        sPtr=baseSection;
    
    /* Process section prefixes of relative path: */
    while(true)
    {
        /* Find next slash in suffix: */
        const char* nextSlashPtr;
        for(nextSlashPtr=pathSuffixPtr;*nextSlashPtr!='\0'&&*nextSlashPtr!='/';++nextSlashPtr)
            ;
        
        /* Navigate section hierarchy: */
        if(nextSlashPtr-pathSuffixPtr==0)
        {
            /* Ignore double slashes */
        }
        else if(nextSlashPtr-pathSuffixPtr==1&&pathSuffixPtr[0]=='.')
        {
            /* Ignore self-reference */
        }
        else if(nextSlashPtr-pathSuffixPtr==2&&pathSuffixPtr[0]=='.'&&pathSuffixPtr[1]=='.')
        {
            /* Go up in the section hierarchy: */
            if(sPtr->parent!=0)
                sPtr=sPtr->parent;
        }
        else
        {
            /* Find subsection name in current section: */
            std::string subsectionName(pathSuffixPtr,nextSlashPtr-pathSuffixPtr);
            Section* ssPtr;
            for(ssPtr=sPtr->firstSubsection;ssPtr!=0&&ssPtr->name!=subsectionName;ssPtr=ssPtr->sibling)
                ;
            
            /* Go down in the section hierarchy: */
            if(ssPtr==0)
            {
                /* Add new subsection to current section: */
                sPtr=sPtr->addSubsection(subsectionName);
                isEdited=true;
            }
            else
                sPtr=ssPtr;
        }
        
        if(*nextSlashPtr=='\0')
            break;
        pathSuffixPtr=nextSlashPtr+1;
    }
    
    /* At this point, sPtr points to the correct section, and pathSuffixPtr is empty: */
    return sPtr;
}

const std::string& ConfigurationFile::retrieveTagValue(const ConfigurationFile::Section* baseSection,const char* relativeTagPath) const
{
    /* Determine the section to start the tag search from: */
    const Section* sPtr;
    const char* tagPathSuffixPtr=relativeTagPath;
    if(tagPathSuffixPtr[0]=='/')
    {
        sPtr=&rootSection;
        ++tagPathSuffixPtr;
    }
    else
        sPtr=baseSection;
    
    /* Process section prefixes of relative tag path: */
    while(true)
    {
        /* Find next slash in suffix: */
        const char* nextSlashPtr;
        for(nextSlashPtr=tagPathSuffixPtr;*nextSlashPtr!='\0'&&*nextSlashPtr!='/';++nextSlashPtr)
            ;
        if(*nextSlashPtr=='\0')
            break;
        
        /* Navigate section hierarchy: */
        if(nextSlashPtr-tagPathSuffixPtr==0)
        {
            /* Ignore double slashes */
        }
        else if(nextSlashPtr-tagPathSuffixPtr==1&&tagPathSuffixPtr[0]=='.')
        {
            /* Ignore self-reference */
        }
        else if(nextSlashPtr-tagPathSuffixPtr==2&&tagPathSuffixPtr[0]=='.'&&tagPathSuffixPtr[1]=='.')
        {
            /* Go up in the section hierarchy: */
            if(sPtr->parent!=0)
                sPtr=sPtr->parent;
        }
        else
        {
            /* Find subsection name in current section: */
            std::string subsectionName(tagPathSuffixPtr,nextSlashPtr-tagPathSuffixPtr);
            Section* ssPtr;
            for(ssPtr=sPtr->firstSubsection;ssPtr!=0&&ssPtr->name!=subsectionName;ssPtr=ssPtr->sibling)
                ;
            
            /* Go down in the section hierarchy: */
            if(ssPtr==0)
                throw TagNotFoundError(relativeTagPath);
            sPtr=ssPtr;
        }
        
        tagPathSuffixPtr=nextSlashPtr+1;
    }
    
    /* At this point, sPtr points to the correct section, and tagPathSuffixPtr is a slash-free tag name: */
    std::string tagName(tagPathSuffixPtr);
    std::list<TagValue>::const_iterator tvIt;
    for(tvIt=sPtr->values.begin();tvIt!=sPtr->values.end()&&tvIt->tag!=tagName;++tvIt)
        ;
    
    /* Return tag value: */
    if(tvIt==sPtr->values.end())
        throw TagNotFoundError(relativeTagPath);
    return tvIt->value;
}

std::string ConfigurationFile::retrieveTagValue(const ConfigurationFile::Section* baseSection,const char* relativeTagPath,const std::string& defaultValue) const
{
    /* Determine the section to start the tag search from: */
    const Section* sPtr;
    const char* tagPathSuffixPtr=relativeTagPath;
    if(tagPathSuffixPtr[0]=='/')
    {
        sPtr=&rootSection;
        ++tagPathSuffixPtr;
    }
    else
        sPtr=baseSection;
    
    /* Process section prefixes of relative tag path: */
    while(true)
    {
        /* Find next slash in suffix: */
        const char* nextSlashPtr;
        for(nextSlashPtr=tagPathSuffixPtr;*nextSlashPtr!='\0'&&*nextSlashPtr!='/';++nextSlashPtr)
            ;
        if(*nextSlashPtr=='\0')
            break;
        
        /* Navigate section hierarchy: */
        if(nextSlashPtr-tagPathSuffixPtr==0)
        {
            /* Ignore double slashes */
        }
        else if(nextSlashPtr-tagPathSuffixPtr==1&&tagPathSuffixPtr[0]=='.')
        {
            /* Ignore self-reference */
        }
        else if(nextSlashPtr-tagPathSuffixPtr==2&&tagPathSuffixPtr[0]=='.'&&tagPathSuffixPtr[1]=='.')
        {
            /* Go up in the section hierarchy: */
            if(sPtr->parent!=0)
                sPtr=sPtr->parent;
        }
        else
        {
            /* Find subsection name in current section: */
            std::string subsectionName(tagPathSuffixPtr,nextSlashPtr-tagPathSuffixPtr);
            Section* ssPtr;
            for(ssPtr=sPtr->firstSubsection;ssPtr!=0&&ssPtr->name!=subsectionName;ssPtr=ssPtr->sibling)
                ;
            
            /* Go down in the section hierarchy: */
            if(ssPtr==0)
                return defaultValue; // Bail out and return default
            else
                sPtr=ssPtr;
        }
        
        tagPathSuffixPtr=nextSlashPtr+1;
    }
    
    /* At this point, sPtr points to the correct section, and tagPathSuffixPtr is a slash-free tag name: */
    std::string tagName(tagPathSuffixPtr);
    std::list<TagValue>::const_iterator tvIt;
    for(tvIt=sPtr->values.begin();tvIt!=sPtr->values.end()&&tvIt->tag!=tagName;++tvIt)
        ;
    
    /* Return tag value: */
    if(tvIt==sPtr->values.end())
        return defaultValue;
    else
        return tvIt->value;
}

std::string& ConfigurationFile::retrieveTagValue(ConfigurationFile::Section* baseSection,const char* relativeTagPath,const std::string& defaultValue)
{
    /* Determine the section to start the tag search from: */
    Section* sPtr;
    const char* tagPathSuffixPtr=relativeTagPath;
    if(tagPathSuffixPtr[0]=='/')
    {
        sPtr=&rootSection;
        ++tagPathSuffixPtr;
    }
    else
        sPtr=baseSection;
    
    /* Process section prefixes of relative tag path: */
    while(true)
    {
        /* Find next slash in suffix: */
        const char* nextSlashPtr;
        for(nextSlashPtr=tagPathSuffixPtr;*nextSlashPtr!='\0'&&*nextSlashPtr!='/';++nextSlashPtr)
            ;
        if(*nextSlashPtr=='\0')
            break;
        
        /* Navigate section hierarchy: */
        if(nextSlashPtr-tagPathSuffixPtr==0)
        {
            /* Ignore double slashes */
        }
        else if(nextSlashPtr-tagPathSuffixPtr==1&&tagPathSuffixPtr[0]=='.')
        {
            /* Ignore self-reference */
        }
        else if(nextSlashPtr-tagPathSuffixPtr==2&&tagPathSuffixPtr[0]=='.'&&tagPathSuffixPtr[1]=='.')
        {
            /* Go up in the section hierarchy: */
            if(sPtr->parent!=0)
                sPtr=sPtr->parent;
        }
        else
        {
            /* Find subsection name in current section: */
            std::string subsectionName(tagPathSuffixPtr,nextSlashPtr-tagPathSuffixPtr);
            Section* ssPtr;
            for(ssPtr=sPtr->firstSubsection;ssPtr!=0&&ssPtr->name!=subsectionName;ssPtr=ssPtr->sibling)
                ;
            
            /* Go down in the section hierarchy: */
            if(ssPtr==0)
            {
                /* Add new subsection to current section: */
                sPtr=sPtr->addSubsection(subsectionName);
                isEdited=true;
            }
            else
                sPtr=ssPtr;
        }
        
        tagPathSuffixPtr=nextSlashPtr+1;
    }
    
    /* At this point, sPtr points to the correct section, and tagPathSuffixPtr is a slash-free tag name: */
    std::string tagName(tagPathSuffixPtr);
    std::list<TagValue>::iterator tvIt;
    for(tvIt=sPtr->values.begin();tvIt!=sPtr->values.end()&&tvIt->tag!=tagName;++tvIt)
        ;
    
    /* Return tag value: */
    if(tvIt==sPtr->values.end())
    {
        /* Add a new tag/value pair: */
        sPtr->values.push_back(TagValue(tagName,defaultValue));
        isEdited=true;
        return sPtr->values.back().value;
    }
    else
        return tvIt->value;
}

ConfigurationFile::ConfigurationFile(const char* sFileName)
    :fileName(sFileName),rootSection(0,""),currentSection(&rootSection),isEdited(false)
{
    /* Load the configuration file: */
    load();
}

ConfigurationFile::~ConfigurationFile(void)
{
}

void ConfigurationFile::load(void)
{
    /* Delete current configuration: */
    while(rootSection.firstSubsection!=0)
    {
        Section* next=rootSection.firstSubsection->sibling;
        delete rootSection.firstSubsection;
        rootSection.firstSubsection=next;
    }
    rootSection.values.clear();
    currentSection=&rootSection;
    
    /* Try opening configuration file: */
    FILE* file=fopen(fileName.c_str(),"rt");
    if(file==0)
        throw FileOpenError(fileName);
    
    /* Read existing configuration file: */
    char line[1024];
    int lineNumber=0;
    while(true)
    {
        /* Read next line from configuration file: */
        fgets(line,1024,file);
        if(feof(file))
            break;
        ++lineNumber;

        /* Remove newline character from end of line: */
        char* lineEndPtr;
        for(lineEndPtr=line;*lineEndPtr!='\n'&&*lineEndPtr!='\0';++lineEndPtr)
            ;
        if(*lineEndPtr=='\0') // Did line fit into read buffer?
            throw MalformedConfigFileError("Line too long",lineNumber,fileName);
        *lineEndPtr='\0';

        /* Skip initial whitespace: */
        char* linePtr;
        for(linePtr=line;*linePtr!='\0'&&isspace(*linePtr);++linePtr)
            ;

        if(*linePtr=='\0'||*linePtr=='#')
            continue; // Ignore empty lines or lines starting with a '#' after whitespace

        /* Read first token from line: */
        char* tokenEndPtr;
        for(tokenEndPtr=linePtr;tokenEndPtr!='\0'&&!isspace(*tokenEndPtr);++tokenEndPtr)
            ;
        *tokenEndPtr='\0';
        if(strcasecmp(linePtr,"section")==0)
        {
            /* Line starts a new section: */
            char* sectionNamePtr;
            for(sectionNamePtr=tokenEndPtr+1;*sectionNamePtr!='\0'&&isspace(*sectionNamePtr);++sectionNamePtr)
                ;
            if(*sectionNamePtr!='\0')
            {
                char* sectionNameEndPtr;
                for(sectionNameEndPtr=sectionNamePtr+1;*sectionNameEndPtr!='\0'&&!isspace(*sectionNameEndPtr);++sectionNameEndPtr)
                    ;
                *sectionNameEndPtr='\0';
                currentSection=currentSection->addSubsection(sectionNamePtr);
            }
            else
                throw MalformedConfigFileError("Missing section name",lineNumber,fileName);
        }
        else if(strcasecmp(linePtr,"endsection")==0)
        {
            /* Line ends a section: */
            if(currentSection->parent!=0)
                currentSection=currentSection->parent;
            else
                throw MalformedConfigFileError("Extra endsection command",lineNumber,fileName);
        }
        else
        {
            /* Line is tag/value pair: */
            char* valuePtr;
            for(valuePtr=tokenEndPtr+1;*valuePtr!='\0'&&isspace(*valuePtr);++valuePtr)
                ;
            currentSection->values.push_back(TagValue(linePtr,valuePtr));
        }
    }

    /* Clean up: */
    fclose(file);
    currentSection=&rootSection;
    isEdited=false;
}

void ConfigurationFile::save(void)
{
    /* Check if the configuration was edited (don't save otherwise): */
    if(isEdited)
    {
        FILE* file=fopen(fileName.c_str(),"wt");
        if(file==0)
            throw FileOpenError(fileName);
        rootSection.save(file,0);
        fclose(file);
        isEdited=false;
    }
}

std::string ConfigurationFile::getCurrentSection(void) const
{
    /* Assemble section prefix string on the way up from the current section: */
    std::string result;
    for(Section* sPtr=currentSection;sPtr!=0;sPtr=sPtr->parent)
        result=sPtr->name+std::string("/")+result;
    return result;
}

void ConfigurationFile::setCurrentSection(const char* newRelativePath)
{
    currentSection=getSection(currentSection,newRelativePath);
}

void ConfigurationFile::list(void) const
{
    /* List all subsections of current section: */
    for(Section* sPtr=currentSection->firstSubsection;sPtr!=0;sPtr=sPtr->sibling)
        printf("%s/\n",sPtr->name.c_str());
    
    /* List all tags in current section: */
    for(std::list<TagValue>::const_iterator tvIt=currentSection->values.begin();tvIt!=currentSection->values.end();++tvIt)
        printf("%s\n",tvIt->tag.c_str());
}
