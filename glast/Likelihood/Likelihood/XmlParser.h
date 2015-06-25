/**
 * @file XmlParser.h
 * @brief Singleton wrapper for xmlBase::XmlParser instance
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/Likelihood/Likelihood/XmlParser.h,v 1.1.1.1 2012/01/17 17:18:14 elwinter Exp $
 */

#ifndef Likelihood_XmlParser_h
#define Likelihood_XmlParser_h

#include "xmlBase/XmlParser.h"

#include "st_facilities/libStApiExports.h"

namespace Likelihood {

/**
 * @class XmlParser
 * @brief Provides a Singleton wrapper a the static instance
 * of xmlBase::XmlParser.
 *
 * @author J. Chiang
 */

class SCIENCETOOLS_API XmlParser {

public:

   /// @return Pointer to the single xmlBase::XmlParser instance.
   static xmlBase::XmlParser * instance() {
      if (s_instance == 0) {
         s_instance = new xmlBase::XmlParser();
      }
      return s_instance;
   }

   static void delete_instance() {
      delete s_instance;
      s_instance = 0;
   }

private:

   static xmlBase::XmlParser * s_instance;

};

// Opaque wrapper since linkage of exported symbols from windows dlls is
// all fouled up.
xmlBase::XmlParser * XmlParser_instance();

} // namespace Likelihood

#endif // Likelihood_XmlParser_h
