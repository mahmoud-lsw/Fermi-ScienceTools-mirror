/**
 * @file SkyConeCut.h
 * @brief Describe an acceptance cone on the sky.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/dataSubselector/dataSubselector/SkyConeCut.h,v 1.4 2011/04/26 04:28:35 jchiang Exp $
 */

#ifndef dataSubselector_SkyConeCut_h
#define dataSubselector_SkyConeCut_h

#include <vector>

#include "astro/SkyDir.h"

#include "dataSubselector/CutBase.h"

namespace dataSubselector {

/**
 * @class SkyConeCut
 * @brief Acceptance cone on the sky.
 */

class SkyConeCut : public CutBase {

public:

   SkyConeCut(double ra, double dec, double radius) 
      : CutBase("SkyCone"), m_ra(ra), m_dec(dec),
        m_coneCenter(astro::SkyDir(ra, dec)), m_radius(radius) {}

   SkyConeCut(const astro::SkyDir & dir, double radius) : 
      CutBase("SkyCone"), m_coneCenter(dir), m_radius(radius) {
      m_ra = m_coneCenter.ra();
      m_dec = m_coneCenter.dec();
   }

   SkyConeCut(const std::string & type, const std::string & unit, 
              const std::string & value);

   virtual ~SkyConeCut() {}

   virtual bool accept(tip::ConstTableRecord & row) const;

   virtual bool accept(const std::map<std::string, double> & params) const;

   virtual CutBase * clone() const {return new SkyConeCut(*this);}

   virtual bool supercedes(const CutBase &) const;
      
   /// @return The corresponding filter string, using the haversine
   /// formula, that can be passed to tip, following the cfitsio
   /// extended filename syntax.  The implementation is from Tom
   /// Stephens' original cutParameters.cxx code.
   virtual std::string filterString() const;

   /// @brief The RA of the cone center (J2000 degrees)
   double ra() const {return m_ra;}

   /// @brief The Dec of the cone center (J2000 degrees)
   double dec() const {return m_dec;}

   /// @brief The cone opening half-angle (degrees)
   double radius() const {return m_radius;}

protected:

   virtual bool equals(const CutBase & rhs) const;

   virtual void getKeyValues(std::string & type, std::string & unit,
                             std::string & value, std::string & ref) const;

private:

   double m_ra;
   double m_dec;
   astro::SkyDir m_coneCenter;
   double m_radius;

   bool accept(double ra, double dec) const;
   void getArgs(const std::string & value, 
                std::vector<std::string> & args) const;

};

} // namespace dataSubselector

#endif // dataSubselector_SkyConeCut_h
