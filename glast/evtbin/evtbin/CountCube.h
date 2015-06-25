/** \file CountCube.h
    \brief Encapsulation of a count map, with methods to read/write using tip.
    \author James Peachey, HEASARC
*/
#ifndef evtbin_CountCube_h
#define evtbin_CountCube_h

#include <string>

#include "evtbin/DataProduct.h"
#include "evtbin/Hist3D.h"

namespace astro {
  class SkyProj;
}

namespace evtbin {

  /** \class CountCube
      \brief Encapsulation of a count map, with methods to read/write using tip.
  */
  class CountCube : public DataProduct {
    public:
      /** \brief Create the count map object.
      */
      CountCube(const std::string & event_file, const std::string & event_table, const std::string & sc_file,
        const std::string & sc_table, double ref_ra, double ref_dec, const std::string & proj,
        unsigned long num_x_pix, unsigned long num_y_pix, double pix_scale, double axis_rot, bool use_lb,
        const std::string & ra_field, const std::string & dec_field, const Binner & energy_binner, const Binner & ebounds, const Gti & gti);

      virtual ~CountCube() throw();

      /** \brief Bin input from input file/files passed to the constructor.
      */
      virtual void binInput();

      /** \brief Bin input from tip table.
          \param begin Table iterator pointing to the first record to be binned.
          \param end Table iterator pointing to one past the last record to be binned.
      */
      virtual void binInput(tip::Table::ConstIterator begin, tip::Table::ConstIterator end);

      /** \brief Write count map file.
          \param creator The value to write for the "CREATOR" keyword.
          \param out_file The output file name.
      */
      virtual void writeOutput(const std::string & creator, const std::string & out_file) const;

    private:
      Hist3D m_hist;
      std::string m_proj_name;
      double m_crpix[2];
      double m_crval[2];
      double m_cdelt[2];
      double m_axis_rot;
      astro::SkyProj * m_proj;
      bool m_use_lb;
      Binner * m_ebounds;
  };

}

#endif
