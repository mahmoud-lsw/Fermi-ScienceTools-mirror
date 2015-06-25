/*
$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/astro/src/IGRField.cxx,v 1.9 2015/01/12 22:54:07 jchiang Exp $
*/
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "igrf_sub/igrf_sub.h"
#include "astro/IGRField.h"

namespace {
   // static instance like singleton 
   astro::IGRField*  model(0);
}

namespace astro {

IGRField& IGRField::Model() { 
   if (model==0) {
      model = new IGRField();
   }
   return *model; 
}

IGRField::IGRField() {
    IGRFf2c::initize_(); 
}

void IGRField::setYear(const float year){
   if(fabs(year-m_year)<0.001) return;
   m_year=year;

   float max_year(2020.);
   if (max_year - m_year < 0.5 && m_year < max_year) {
      std::cout << "WARNING: Requested year, " << m_year 
                << ", is within six months of the "
                << "end of valid range (1900-" << max_year << ") "
                << "for the IGRF-12 model." << std::endl;
   }
   if (m_year >= max_year) {
      if (!std::getenv("IGNORE_IGRF_BOUNDARY")) {
         std::ostringstream message;
         message << "Requested year, " << m_year 
                 << ", is outside the valid range of 1900-" << max_year;
         throw std::runtime_error(message.str());
      }
      m_year = max_year;
   }
   if (m_year <= 1900.) {
      std::ostringstream message;
      message << "Requested year, " << year 
              << ", is outside the valid range of 1900-2015";
      throw std::runtime_error(message.str());
   }
   IGRFf2c::feldcof_(&m_year,&m_dipoleMoment);
}
      
int IGRField::compute(const float latitude,const float longitude,const float altitude,const float year){

   float earth_radius= 6371.2 ;

   float b0,rr0;
   IGRFf2c::integer icode;
   IGRFf2c::logical value;
   float stps=0.05;
   float bdel=0.001; 

   if(year>0) Model().setYear(year);
   m_latitude=latitude;
   m_longitude=longitude;
   m_altitude=altitude;
   
   IGRFf2c::feldg_(&m_latitude,&m_longitude,&m_altitude,&m_bNorth,&m_bEast,&m_bDown,&m_bAbs);
   IGRFf2c::shellg_(&m_latitude,&m_longitude,&m_altitude,&m_dipoleMoment,&m_L,&icode,&b0) ;
   
   IGRFf2c::findb0_(&stps, &bdel, &value, &m_bEquator, &rr0);
   if(value==0) m_bEquator=m_dipoleMoment/(m_L*m_L*m_L);
   
   m_B=m_bAbs/m_bEquator;   
   
// conversion from dipole moment in IGRF G*R_earth^3 to rigidity in GV:
//                                                R_earth   km->cm  statvolt -> gigavolt
   float rigidity_const= 0.25 * m_dipoleMoment * earth_radius *  1e5 *  300. / 1e9 ;
   
// the invariant latitude and the rigidity cutoff
// as described in Smart and Shea, Adv. Space Res. 36 (2005) p.2012

   m_invariantLat=m_L>=1 ? acos(sqrt(1./m_L)): 0.;
   m_rigidityCutoff=rigidity_const/(m_L*m_L);
   
// and alternatively the r, lambda coordinates 
// from http://www.spenvis.oma.be/spenvis/help/background/magfield/rlambda.html 
// gives identical rigidity cutoffs if rigidity cutoff is calculated as 
// rigidity_const * cos(m_lambda)^4 / m_R^2

// needs an approximation to the solution rl(b) of the equation b^2*rl^6+3*rl-4==0 in the range b=1..10
   double Bdoub = (double)m_B;
   double rl= pow(Bdoub,-0.215108)*(1.-0.020551*(Bdoub-1.)+0.0008148*(Bdoub-1.)*(Bdoub-1.));
   m_R= rl*m_L;
   m_lambda= (rl<=1) ? acos(sqrt(rl)) : 0.;

   return icode;

}

const IGRF_data & IGRField::igrf_data(const std::string & filename) {
   std::map<std::string, astro::IGRF_data>::const_iterator it 
      = m_IGRF_data_map.find(filename);
   if (it != m_IGRF_data_map.end()) {
      return it->second;
   }
   m_IGRF_data_map.insert(std::make_pair(filename, IGRF_data(filename)));
   return m_IGRF_data_map[filename];
}


}
