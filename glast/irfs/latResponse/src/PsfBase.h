/**
 * @file PsfBase.h
 * @brief PsfBase class declaration.
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/irfs/latResponse/src/PsfBase.h,v 1.1.1.3.2.4 2015/04/26 06:53:49 jasercio Exp $
 */

#ifndef latResponse_PsfBase_h
#define latResponse_PsfBase_h

#include <string>

#include "irfInterface/IPsf.h"

namespace latResponse {

/**
 * @class PsfBase
 *
 */

class PsfBase : public irfInterface::IPsf {

public:

   PsfBase(const std::string & fitsfile, bool isFront,
           const std::string & extname="RPSF", 
           const std::string & scaling_extname="PSF_SCALING_PARAMS");

   PsfBase(const PsfBase & rhs);

   PsfBase & operator=(const PsfBase &);

   virtual ~PsfBase() {}

   virtual double scaleFactor(double energy, bool isFront) const;

   typedef std::vector<irfInterface::AcceptanceCone *> AcceptanceConeVector_t;

protected:

   virtual double scaleFactor(double energy) const;

private:

   // PSF scaling parameters
   double m_par0;
   double m_par1;
   double m_index;

   // store all of the PSF parameters
   std::vector<double> m_psf_pars;

   void readScaling(const std::string & fitsfile, bool isFront,
                    const std::string & extname);
};

} // namespace latResponse

#endif // latResponse_PsfBase_h
