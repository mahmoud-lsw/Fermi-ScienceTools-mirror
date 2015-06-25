/**
 * @file Ft2File.h
 * @brief Declaration of FT1 file abstraction
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/fitsGen/fitsGen/Ft2File.h,v 1.3 2008/03/13 22:57:56 jchiang Exp $
 */

#ifndef fitsGen_Ft2File_h
#define fitsGen_Ft2File_h

#include "fitsGen/FtFileBase.h"

namespace fitsGen {

/**
 * @class Ft2File
 * @brief Abstraction/interface layer for using tip to write FT2
 * files.
 *
 * @author J. Chiang
 */

class Ft2File : public FtFileBase {

public:

   Ft2File(const std::string & outfile,
           long nrows=0, 
           const std::string & extname="SC_DATA",
           const std::string & templateFile="ft2.tpl")
      : FtFileBase(outfile, nrows) {
      init(templateFile, extname);
   }

   /// Fill the RA_SCZ, DEC_SCZ, RA_SCX, DEC_SCX values and the
   /// quaternions, QSJ_1, QSJ_2, QSJ_3, QSJ_4
   void setScAxes(double ra_scz, double dec_scz,
                  double ra_scx, double dec_scx);

};

} // namespace fitsGen

#endif // fitsGen_Ft2File_h
