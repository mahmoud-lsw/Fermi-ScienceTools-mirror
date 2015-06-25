#ifndef FluxExecption_H_
#define FluxExecption_H_

/** 
* \class FluxException
*
* \brief Define FATAL MACRO, which throws exception with error message
* \author Toby Burnett
* 
* $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/flux/FluxException.h,v 1.1.1.1 2003/07/29 18:22:14 burnett Exp $
*/

#include <sstream>


#define FATAL_MACRO(output)\
    do{std::stringstream message; \
    message <<__FILE__<<":"<<__LINE__<<": "<<output<<'\0';\
    std::cerr << message.str() << std::endl; \
    throw(message.str());}while(0)

#define WARNING_MACRO(output)\
std::cerr << output << std::endl;

#endif // _ERROR_H_

