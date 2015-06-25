/**
 * @file Primack05.h
 * @brief Encapsulation of Primack, Bullock, Somerville (2005) astro-ph 0502177
 * EBL model.
 * 
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/celestialSources/eblAtten/src/Primack05.h,v 1.1.1.2 2011/03/20 19:24:56 elwinter Exp $
 */

#include <vector>

namespace IRB {

/**
 * @class Primack05
 *
 */

class Primack05 {

public:

   static Primack05 & instance();
   
   float value(float energy, float redshift) const;

protected:

   static Primack05 * s_instance;

   Primack05();

private:

   std::vector<float> m_zvalue;
   std::vector<float> m_evalue;
   std::vector< std::vector<float> > m_tauTables;

   void fillTauTables();

};

} // namespace IRB
