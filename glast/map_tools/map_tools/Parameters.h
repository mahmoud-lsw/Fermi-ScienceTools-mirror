/**
* @file Parameters.h
* @brief Tool Input Parameter Reader base class
*
* $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/map_tools/map_tools/Parameters.h,v 1.15 2006/02/04 17:47:54 jchiang Exp $
*/
#ifndef MAP_TOOLS_PARAMETERS_H
#define MAP_TOOLS_PARAMETERS_H 

#include <string>

#include "hoops/hoops_prompt_group.h"

namespace map_tools {
/**
* @class Parameters
* @brief Input reader base class for tools
*
* It uses hoops to read parameters from the par file.
* The description of pil format is available at
* <a href="http://www-glast.slac.stanford.edu/sciencetools/userInterface/doc/pil.pdf">PIL user
* manual</a>.
*
* @author Toby Burnett [originally from Sandhia Bansall]
*
*/
class Parameters 
{
public:
    //! @brief ctor for command line args
    Parameters( int argc, char *argv[]);

    //! @brief ctor with hoops already set up
    Parameters( hoops::ParPromptGroup& par);

    ~Parameters() throw() ;

    template< typename T>
        T getValue(const std::string & name) const { return m_par[name];}

     template <typename T>
   T getValue(const std::string & name, const T & default_value) const {
      try {
         return getValue<T>(name);
      } catch ( ...) {
         return default_value;
      }
   }
   //! @brief bracket operator for numeric values only
   double operator[](const std::string& name)const {return m_par[name];} 

    // Accessor Methods
    const std::string &inputFile() const   { return m_inFile; }
    const std::string &table_name() const   { return m_table_name; }
    const std::string &filter() const      { return m_filter; }
    const std::string &outputFile() const  { return m_outFile; }
    bool verboseMode()  const            { return m_verboseMode; }
    bool clobber()      const            { return m_clobber; }
    short chatter()     const            { return m_chatter; }
private:
    void setup();
     hoops::ParPromptGroup & m_par;
    // Data Members
    std::string   m_inFile;
    std::string   m_outFile;
    std::string   m_filter;
    std::string   m_table_name;
    bool        m_verboseMode;
    bool        m_clobber;
    int         m_chatter;
    bool      m_own_ppg;
    bool    m_debug;
};
} // namespace map_tools
#endif
