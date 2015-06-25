/** 
* @file SpectrumFactoryLoader.h
* @brief Load the external spectrum factory objects
*
*  $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/celestialSources/SpectrumFactoryLoader.h,v 1.2 2005/04/27 22:11:40 burnett Exp $
*/

#ifndef SpectrumFactoryLoader_h
#define SpectrumFactoryLoader_h

#include <vector>
#include <string>

class ISpectrumFactory;

class SpectrumFactoryLoader {
public:
    /// ctor does the work
    SpectrumFactoryLoader();
    /// access to a list of the names that were loaded
    std::vector<std::string> names()const{return m_names;}
private:
    void load(ISpectrumFactory&);
    std::vector<std::string> m_names;
};

#endif
