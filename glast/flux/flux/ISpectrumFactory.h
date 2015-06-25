/** 
* @file ISpectrumFactory.h
* @brief declare ISpectrumFactory interface
*
*  $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/flux/ISpectrumFactory.h,v 1.2 2005/04/26 17:19:35 burnett Exp $
*/

#ifndef flux_ISpectrumFactory_h
#define flux_ISpectrumFactory_h

class ISpectrum;

/** 
* \class ISpectrumFactory
*
* \brief This is an abstract base class for the SpectrumFactory template classes.
* 
*/
class ISpectrumFactory
{
public:
    
    /// the only thing this does: make an associated Spectrum object
    virtual ISpectrum* instantiate(const std::string& params)const=0;
    
    /// a dummy to force creation of the spectrum object.
    virtual void addRef()const=0;

    //! access to the name
    virtual std::string name()const=0;
    
};


#endif 
