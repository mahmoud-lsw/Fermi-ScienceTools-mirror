/**
 * @file Util.h
 * @brief Utility functions used by response function classes and others.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/irfs/irfUtil/irfUtil/Util.h,v 1.10 2015/01/07 22:55:21 jchiang Exp $
 */

#ifndef irfUtil_Util_h
#define irfUtil_Util_h

#include <string>
#include <map>
#include <vector>

namespace irfUtil {

/**
 * @class Util
 *
 */

class Util {

public:

   virtual ~Util() {}

   static void getCaldbFile(const std::string & detName, 
                            const std::string & respName,
                            const std::string & version,
                            std::string & filename,
                            long & extnum,
                            const std::string & telescope="GLAST",
                            const std::string & instrument="LAT",
                            const std::string & filter="NONE",
                            const std::string & date="-",
                            const std::string & time="-");

   static void get_event_class_mapping(std::map<std::string, unsigned int> 
                                       & event_class_mapping);

   static void get_event_type_mapping(const std::string & event_class,
                                      std::map<std::string,
                                      std::pair<unsigned int, std::string> >
                                      & event_type_mapping); 

   static void get_event_type_mapping(const std::string & event_class,
                                      std::map<std::string,
                                      std::pair<unsigned int, std::string> >
                                      & event_type_mapping,
                                      std::vector<std::string> & partitions,
                                      std::map<std::string, unsigned int> 
                                      & bitmask_by_partition);

protected:

   Util() {}

};

} // namespace irfUtil

#endif // irfUtil_Util_h
