/** @file FluxMgr.h
    @brief declaration of FluxMgr

 $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/flux/FluxMgr.h,v 1.19 2007/08/15 18:46:42 burnett Exp $

  */
#ifndef FLUX_MGR_H
#define FLUX_MGR_H

/** 
* \class FluxMgr
*
* \brief The point of entry for interfacing with the flux package.
* holds methods for creating sources and sending new particles, 
* and methods for interfacing with the satellite position, and 
* setting the position variables. It is instantiated with
* the names of the xml files to be used as input to the xml parser.
* 
* 
*/

#include "astro/GPS.h"

#include "FluxSource.h"

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/XmlParser.h"
#include "ISpectrumFactory.h"
#include <map>
#include <list>
#include <string>

using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;

class FluxMgr 
{
    
public:
    
    /// ctor for multiple XML documents
    FluxMgr(const std::vector<std::string>& fileList, std::string dtd="");
    
    ~FluxMgr();
    
    /// create and return a source by name.
    EventSource* source(std::string name);

    /// create a composite source from the list of names
    EventSource* compositeSource(std::vector<std::string> names);

    
    /// access to the source list
    std::list<std::string> sourceList() const;
    
    /// set the target area
    void setArea(double area);
    
    /// generate some test output
    void test(std::ostream& out, std::string source_name, int count);

    /// set the desired pointing history file to use:
    void setPointingHistoryFile(std::string fileName);

    ///this should return the source file names, along with the contained sources.
    std::vector<std::pair< std::string ,std::list<std::string> > > sourceOriginList() const;
    
    void addFactory(std::string name, const ISpectrumFactory* factory );
    
    /// set the expansion factor for the orbit (-1) = random
    void setExpansion (double p);

    /// pass a specific amount of time
    void pass ( double t);

    /// Get the time as held by GPS
    double time () const;

    /// Set the time
    void setTime(double newtime);

    /// synch satellite location with current time
    void synch ();
    
    /// set the sample interval
    void sampleintvl ( /*GPStime*/double t );
    
    /// get the current satellite location
    std::pair<double,double> location();
    
    CLHEP::HepRotation transformToGlast(double seconds, astro::GPS::CoordSystem index);

    ///this sets the rocking mode in GPS.
    std::vector<double> setRockType(astro::GPS::RockType rockType, double rockAngle);
    
    /// Set an alignment rotation:
    /// @param qx,qy, qz rotation angles (degrees) about coordinate axes -- assume small, < 1 deg
    /// @param misalign if true, apply to incoming; if false, apply as correction to coordinate transformation
    ///
    void setAlignmentRotation(double qx, double qy, double qz, bool misalign);

    /// set an offset for generating source id numbers, return previous value
    int setIdOffset(int id);

    /** set a cone to filter incoming (galactic) data
    @param ra,dec center of cone, equatorial coords in degrees
    @param radius radius of cone, degrees
    */
    void setFilterCone(double ra, double dec, double radius);

private:
    
    /// source library lookup.  Each source is uniquely identified
    /// by its "name" attribute because "name" is of type ID
    DOMElement* getLibrarySource(const std::string& id);
    
    
    void defaultFile();
    void init(const std::vector<std::string>& fileList);
    
    EventSource* 
    getSourceFromXML(const DOMElement* src);
    
    DOMDocument* m_library_doc;
    
    DOMElement*     s_library;
    
    std::vector<DOMDocument*> m_library_doclist;
    
    std::vector<DOMElement*>     s_librarylist;
    
    /// list of sources for easy lookup
    //std::map<std::string, DOM_Element > m_sources;
    std::map<std::string, std::pair<DOMElement* ,std::string> > m_sources;

    /// internal routine that creates the document
    std::string  writeXmlFile( const std::vector<std::string>& fileList);
    
    /// filename for dtd
    std::string m_dtd;
};
#endif
