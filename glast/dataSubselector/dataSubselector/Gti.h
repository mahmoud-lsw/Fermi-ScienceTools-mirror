/**
 * @file Gti.h
 * @brief Handle GTIs. Derives from evtbin::Gti, adding a constructor and 
 * an accept() method.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/dataSubselector/dataSubselector/Gti.h,v 1.8 2011/04/26 04:28:35 jchiang Exp $
 */

#ifndef dataSubselector_Gti_h
#define dataSubselector_Gti_h

#include "evtbin/Gti.h"

namespace tip {
   class Table;
}

namespace dataSubselector {

/**
 * @class Gti
 *
 */

class Gti : public evtbin::Gti {

public:

   Gti() : evtbin::Gti() {}

   Gti(const std::string & filename, const std::string & extension="GTI") 
      : evtbin::Gti(filename, extension) {}

   Gti(const tip::Table & gtiTable);

   Gti(const evtbin::Gti & gti);

   bool accept(double time) const;

   bool accept2(double time) const;

   void writeExtension(const std::string & filename) const;

   /// @brief Given a start and stop time this method recomputes the
   /// the GTIs.
   /// @return A new Gti object with the new good-time intervals.
   /// @param start The start time of the time range cut (MET seconds)
   /// @param stop The stop time of the time range cut (MET seconds)
   Gti applyTimeRangeCut(double start, double stop) const;

   /// @return The minimum lower bound of the GTIs (MET seconds)
   double minValue() const;

   /// @return The maximum upper bound of the GTIs (MET seconds)
   double maxValue() const;

};

} // namespace dataSubselector

#endif // dataSubselector_Gti_h
