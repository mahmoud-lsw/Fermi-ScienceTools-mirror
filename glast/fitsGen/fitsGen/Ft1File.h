/**
 * @file Ft1File.h
 * @brief Declaration of FT1 file abstraction
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/fitsGen/fitsGen/Ft1File.h,v 1.9 2005/12/23 19:51:50 jchiang Exp $
 */

#ifndef fitsGen_Ft1File_h
#define fitsGen_Ft1File_h

#include "astro/JulianDate.h"

#include "fitsGen/FtFileBase.h"

namespace fitsGen {

/**
 * @class Ft1File
 * @brief Abstraction/interface layer for using tip to write FT1
 * files.
 *
 * @author J. Chiang
 */

class Ft1File : public FtFileBase {

public:

   Ft1File(const std::string & outfile,
           long nrows=0,
           const std::string & extname="EVENTS",
           const std::string & templateFile="ft1.tpl");

   virtual ~Ft1File();

   virtual void close();

private:

   void verifyObsTimes();

};

} // namespace fitsGen

#endif // fitsGen_Ft1File_h
