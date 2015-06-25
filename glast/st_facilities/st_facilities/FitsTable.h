/**
 * @file FitsTable.h
 * @brief Class to manage tabular IRF data which are assumed to be a
 * function of log10(E/MeV) and cos(theta) and are read from a FITS
 * file.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/st_facilities/st_facilities/FitsTable.h,v 1.1 2012/12/10 10:17:45 cohen Exp $
 */

#ifndef st_facilities_FitsTable_h
#define st_facilities_FitsTable_h

#include <map>
#include <vector>

namespace tip {
   class Table;
}

namespace st_facilities {

class Bilinear;

class FitsTable {

public:

   FitsTable(const std::string & filename,
             const std::string & extname,
             const std::string & tablename,
             size_t nrow=0);

   FitsTable();

   FitsTable(const FitsTable & other);

   ~FitsTable();
      
   /// @brief lookup a value from the table
   /// @param logenergy log10(energy)
   /// @param costh cos(theta)
   /// @param interpolate [true] if true, make linear
   /// interpolation. Otherwise, return value for given cell.
   double value(double logenergy, double costh, bool interpolate=true) const;
    
   double maximum() const {
      return m_maxValue;
   }
   
   double minCosTheta() const {
      return m_minCosTheta;
   }

   static void getVectorData(const tip::Table * table,
                             const std::string & fieldName,
                             std::vector<double> & values,
                             size_t nrow=0);

   void getValues(std::vector<double> & values) const;

   void getCornerPars(double logE, double costh, 
                      double & tt, double & uu,
                      std::vector<double> & cornerEnergies,
                      std::vector<double> & cornerPars) const;

   const std::vector<double> & logEnergies() const {
      return m_logEnergies;
   }

   const std::vector<double> & costhetas() const {
      return m_mus;
   }

   // bounds in log10(E)
   const std::vector<double> & ebounds() const {
      return m_ebounds;
   }

   // bounds in cos(theta)
   const std::vector<double> & tbounds() const {
      return m_tbounds;
   }

   double getPar(size_t ilogE, size_t icosth) const;

   void setPar(size_t ilogE, size_t icosth, double par);
   
protected:

   /// Disable copy assignment operator.
   FitsTable & operator=(const FitsTable &) {
      return *this;
   }

private:

   Bilinear * m_interpolator;

   std::vector<double> m_logEnergies; 
   std::vector<double> m_mus; 
   std::vector<double> m_values;

   std::vector<double> m_ebounds;
   std::vector<double> m_tbounds;
   
   double m_minCosTheta;

   double m_maxValue;

};

} // namespace st_facilities

#endif // st_facilities_FitsTable_h
