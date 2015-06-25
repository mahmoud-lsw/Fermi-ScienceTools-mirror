/**
 * @file Cuts.h
 * @brief Handle data selections and DSS keywords for packages outside of
 * dataSubselector.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/dataSubselector/dataSubselector/Cuts.h,v 1.47 2015/02/20 00:09:27 jchiang Exp $
 */

#ifndef dataSubselector_Cuts_h
#define dataSubselector_Cuts_h

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "dataSubselector/GtiCut.h"
#include "dataSubselector/RangeCut.h"

namespace tip {
   class ConstTableRecord;
   class Header;
   class Table;
}

namespace dataSubselector {

class BitMaskCut;
class Gti;
class GtiCuts;

/**
 * @class Cuts
 * @brief Handle data selections and Data Subspace (DSS) keywords for
 * packages outside of dataSubselector.
 * @author J. Chiang
 *
 */

class Cuts {

public: 

   Cuts() {}

   /// @brief This constructor reads the data selections from the event 
   ///        extension header.
   /// @param eventFile FITS file name.
   /// @param extension Name of the HDU containing the cuts.
   /// @param check_columns DSS range cuts may be specified to 
   ///        have a functional dependence on other column values,
   ///        e.g., to specify an energy dependent zenith angle cut,
   ///        the DSTYPn keyword value would be
   ///        "ZENITH_ANGLE(ENERGY)".  Currently, these classes only
   ///        support simple range cuts that do not have these
   ///        functional dependences.  If check_columns is true, an
   ///        exception will be thrown if a non-bare column name is
   ///        encountered as the DSS type.  If false, the range cut
   ///        will be stored, but if applied to a tip::TableRecord,
   ///        tip will throw an exception.
   /// @param skipTimeRangeCuts If true, time range cuts in the DSS
   ///        keywords will not be read in.  Since the GTI extension
   ///        should already include any time range cuts, explicit
   ///        time range cuts are not needed.
   Cuts(const std::string & eventFile, const std::string & extension,
        bool check_columns=true, bool skipTimeRangeCuts=false,
        bool skipEventClassCuts=false);

   /// @brief This constructor reads in a vector of eventFiles, verifying
   ///        that the non-GTI cuts are the same in all files, and merging
   ///        the GTIs from the various files into a single Gti object.
   Cuts(const std::vector<std::string> & eventFiles,
        const std::string & extension,
        bool check_columns=true,
        bool skipTimeRangeCuts=false,
        bool skipEventClassCuts=false);

   /// A copy constructor is needed since there are pointer data members.
   Cuts(const Cuts & rhs);

   /// Deletes the CutBase pointers.
   ~Cuts();

   /// Copy assignment operator
   Cuts & operator=(const Cuts & rhs);

#ifndef SWIG
   /// @brief True if the data in the row passes all of the cuts.
   /// @param row A row of FITS binary table.
   bool accept(tip::ConstTableRecord & row) const;

   /// @brief True if the data in the map passes all of the cuts.
   /// @param params A map of keyword-value pairs for use in contexts
   ///        where filtering is desired but the data do not exist in
   ///        FITS format.
   bool accept(const std::map<std::string, double> & params) const;

   /// @brief This method will add a cut if it is not equal to
   ///        and if it is not superceded by an existing cut.
   ///        If the added cut supercedes an existing cut, that cut
   ///        will be replaced.
   /// @return The current number of cuts stored.
   /// @param colname The column name of the value to be selected on.
   /// @param unit The units of this quantity, e.g., 's', 'deg', 'MeV'
   /// @param minVal The minimum value to be accepted by this cut.
   /// @param maxVal The maximum value to be accepted by this cut.
   /// @param type Specifies whether the cut range is open (bounded 
   ///        from below or above) or closed.
   /// @param indx If the column is a vector, indx is the index of the
   ///        array to be cut on.  FITS follows the FORTRAN convention
   ///        that the first element is indx=1.
   unsigned int addRangeCut(const std::string & colname, 
                            const std::string & unit,
                            double minVal, double maxVal,
                            RangeCut::IntervalType type=RangeCut::CLOSED,
                            unsigned int indx=0);

   /// @brief Add a GTI cut.  Here an existing GTI extension is 
   ///        read from a FITS file as a tip::Table.
   /// @return The current number of cuts stored.
   /// @param gtiTable The GTI extension read using tip.
   unsigned int addGtiCut(const tip::Table & gtiTable);

   /// @brief Add a GTI cut using a Gti object.
   /// @return The current number of cuts stored.
   /// @param gti A dataSubselector::Gti object (NB: This class derives
   ///        from evtbin::Gti).
   unsigned int addGtiCut(const Gti & gti);

   /// @brief An acceptance cone on the sky.
   /// @return The current number of cuts stored.
   /// @param ra Right Ascension of cone center in (J2000 degrees).
   /// @param dec Declination of cone center (J2000 degrees).
   /// @param radius Cone half-opening angle (degrees)
   unsigned int addSkyConeCut(double ra, double dec, double radius);

   /// @brief (single) Bit mask cut, usually used to make event
   /// selection based on bit map in the EVENT_CLASS column
   /// @return The current number of cuts stored.
   /// @param colname Name of column to which the bit-mask is applied
   /// @param bitPosition Position of the bit (corresponding to the
   /// desired event class)
   /// @param pass_ver Pass version of IRFs, e.g., "P7V6".
   ///        This should be left blank for pre-Pass 7 IRFs.
   unsigned int addBitMaskCut(const std::string & colname,
                              unsigned int bitPosition,
                              const std::string & pass_ver="");

   /// @brief Not a real cut since it does not apply to any columns
   /// in FT1.  This cut is used to keep track of the IRF version
   /// information so that it is automatically propagated to derived
   /// data products such as filtered FT1 files, counts maps, etc..
   /// @return The current number of cuts stored.
   /// @param colname Name of IRF property that is versioned.
   /// @param version The version string.
   unsigned int addVersionCut(const std::string & colname,
                              const std::string & version);

   unsigned int addCut(const CutBase & newCut) {
      m_cuts.push_back(newCut.clone());
      return m_cuts.size();
   }

   /// @brief Merge RangeCuts, taking care of intersecting intervals, 
   /// open-ended ranges, etc.
   unsigned int mergeRangeCuts();

   /// @brief Remove all RangeCuts with a given column name.
   /// @param colname The name of the column corresponding to the
   ///        RangeCuts to be removed.
   /// @param removedCuts Vector of removed RangeCuts.  The calling routine
   ///        is responsible for deleting these cuts.
   unsigned int removeRangeCuts(const std::string & colname,
                                std::vector<RangeCut *> & removedCuts);

   /// @brief Remove specified VersionCut
   /// @param colname Column name of VersionCut to be removed.
   unsigned int removeVersionCut(const std::string & colname);

   /// @brief Add the DSS keywords to a write-enabled tip::Header.
   /// @param header The reference returned from the tip::Table::getHeader
   ///        method where the table has been opened using the
   ///        tip::IFileSvc::editTable method.
   void writeDssKeywords(tip::Header & header) const;

   /// @brief Add the time-related DSS keywords to a write-enabled
   ///        tip::Header.
   /// @param header The reference returned from the tip::Table::getHeader
   ///        method where the table has been opened using the
   ///        tip::IFileSvc::editTable method.
   void writeDssTimeKeywords(tip::Header & header) const;

   /// @brief Add the first GTI extension in the m_cuts vector to
   ///        the specified FITS file.  It is tacitly assumed that
   ///        there is only one GTI extension.  If a GTI extension
   ///        already exists in the FITS file, then nothing is done.
   /// @param filename The FITS file name.
   /// @todo This method should check if the existing extension is
   /// compatible with the one to be written. If it is, then the old
   /// one is should be replaced.  If it is not, then a runtime
   /// exception should be thrown.
   void writeGtiExtension(const std::string & filename) const;

   /// @brief The number of cuts contained herein.
   unsigned int size() const {return m_cuts.size();}

   /// @brief Do a member-wise comparison of each cut in m_cuts.
   /// The ordering of cuts must be the same.
   bool operator==(const Cuts & rhs) const;

   bool operator!=(const Cuts & rhs) const {
      return !operator==(rhs);
   }

   /// @brief Do a member-wise comparison of each cut, but skip GTIs
   bool compareWithoutGtis(const Cuts & rhs) const;

   /// @brief Write a summary of the cuts to the output stream.  This
   /// summary contains the same information as the DSS keywords.
   /// @param stream Where the output should be directed.
   /// @param suppressGtis Suppress printing of the GTI data if it exists.
   void writeCuts(std::ostream & stream, bool suppressGtis=false) const;

   /// @return A reference to the i-th cut.
   const CutBase & operator[](unsigned int i) const {return *(m_cuts.at(i));}

   std::string filterString() const;

   /// @return A vector of pointers to the GtiCuts that are present.
   void getGtiCuts(std::vector<const GtiCut *> & gtiCuts);
#endif

   /// @return The irfs named in the DSS keywords.
   const std::string & irfName() const;

   /// @return IRFs implied by the bit position and Pass version in
   /// the BIT_MASK DSS keywords, using (as specified by the FSSC),
   /// the most recent version of the corresponding IRFs in
   /// $CALDB/bcf/irf_index.fits.
   std::string CALDB_implied_irfs() const;

   /// @brief Append the event type partition to the irfs_name.  The
   /// event type partition is determined from the EVENT_TYPE
   /// BitMaskCut value. For FB, the result is just the generic irfs
   /// name, for PSF and EDISP, " (PSF)" and " (EDISP)" are appended.
   /// The addition of the event type partition info is needed to
   /// index the desired irfs from irfLoader.
   /// @param irfs_name The generic irfs name, e.g., P8R2_SOURCE_V6.
   void append_event_type_partition(std::string & irfs_name) const;

#ifndef SWIG
   /// This will set the BitMaskCut and m_pass_ver value based 
   /// on the IRF name (assuming it is Pass 7 or later).
   /// If BitMaskCut is already set, this will throw an exception.
   void setIrfs(const std::string & irfName);

   /// @return A new Cuts object. This static function checks that all
   /// of the non-GtiCuts are the same in each element of the input
   /// vector, copies those, then merges the GtiCuts.  The returned
   /// object is then equivalent to the composition of all of the Cuts
   /// in the vector.
   static Cuts mergeGtis(std::vector<Cuts> & cuts_vector);

   /// @return True if the cut is a time range or GTI cut.
   static bool isTimeCut(const CutBase & cut);

   static void checkIrfs(const std::string & infile,
                         const std::string & extname,
                         const std::string & irfs);

   RangeCut * conversionTypeCut() const;

#endif

   BitMaskCut * bitMaskCut(const std::string & colname="EVENT_CLASS") const;

   std::vector<BitMaskCut *> bitMaskCuts() const;

   const std::string & pass_ver() const {
      return m_pass_ver;
   }

   const bool post_P7() const {
      return m_post_P7;
   }

   /// Replace or add an existing BitMaskCut applied to the same column.
   void setBitMaskCut(BitMaskCut * bitMaskCut);

   static void extract_irf_versions(const std::string & irf_name,
                                    std::string & pass_ver,
                                    std::string & irf_ver);

   static void read_bitmask_mapping(std::map<std::string,
                                    unsigned int> & irfs);

private:

   std::vector<CutBase *> m_cuts;

   std::string m_irfName;

   std::string m_pass_ver;
   bool m_post_P7;

   unsigned int parseColname(const std::string & colname,
                             std::string & col) const;

   bool hasCut(const CutBase * newCut) const;

   unsigned int find(const CutBase * cut) const;

   /// @brief Add a cut. The passed cut will not be added if an 
   ///        existing cut supercedes it, but it will be deleted.  If
   ///        added, this cut will be deleted by the destructor ~Cut().
   unsigned int addCut(CutBase * newCut);

   /// @brief Remove all DSS keywords (and NDSKEYS) from
   ///        the header.
   void removeDssKeywords(tip::Header & header) const;

   void read_pass_ver(const std::string & infile, const std::string & ext);

   void set_irfName(const std::string & infile, const std::string & ext);
};

} // namespace dataSubselector

#endif // dataSubselector_Cuts_h
