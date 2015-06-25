/** \file PulsarDb.h
    \brief Interface for PulsarDb class.
    \authors Masaharu Hirayama, GSSC,
             James Peachey, HEASARC/GSSC
*/
#ifndef pulsarDb_PulsarDb_h
#define pulsarDb_PulsarDb_h

#include <list>
#include <map>
#include <stdexcept>
#include <string>
#include <set>

#include "pulsarDb/EphStatus.h"
#include "pulsarDb/OrbitalEph.h"
#include "pulsarDb/PulsarEph.h"

#include "tip/FileSummary.h"
#include "tip/Table.h"
#include "tip/TipFile.h"

namespace tip {
  class Extension;
  class Header;
}

namespace pulsarDb {

  /** \class IEphFactory
      \brief Base class to define the interface of ephemeris factory classes.
             Note: Template parameter, EPHTYPE, must be either PulsarEph or OrbitalEph.
  */
  template <typename EPHTYPE>
  class IEphFactory {
    public:
      /// \brief Destruct this IEphFactory object.
      virtual ~IEphFactory() {}

      /** \brief Create an object of a subclass of PulsarEph or OrbitalEph from a FITS row.
          \param record FITS row from which ephemeris parameters are to be read.
          \param header FITS header to read other information if necessary.
      */
      virtual EPHTYPE * create(const tip::Table::ConstRecord & record, const tip::Header & header) const = 0;
  };

  /** \class EphFactory
      \brief Concrete class to create a PulsarEph object or an OrbitalEph object.
             Note: The first emplate parameter, EPHTYPE, must be either PulsarEph or OrbitalEph.
                   The second template parameter, EPHSTYLE, must be one of their subclasses.
  */
  template <typename EPHTYPE, typename EPHSTYLE>
  class EphFactory: public IEphFactory<EPHTYPE> {
    public:
      /** \brief Create an object of a subclass of PulsarEph or OrbitalEph from a FITS row.
          \param record FITS row from which ephemeris parameters are to be read.
          \param header FITS header to read other information if necessary.
      */
      virtual EPHTYPE * create(const tip::Table::ConstRecord & record, const tip::Header & header) const {
        return new EPHSTYLE(record, header);
      }

      /// \brief Return a singleton factory for a PulsarEph object or a OrbitalEph object.
      static EphFactory & getFactory() {
        static EphFactory s_factory;
        return s_factory;
      }

    private:
      /// \brief Construct a EphFactory object.
      EphFactory() {};
  };

  /** \class PulsarDb
      \brief Abstraction providing access to pulsar ephemerides database.
  */
  class PulsarDb {
    public:
      typedef std::vector<tip::Table *> TableCont;
      typedef std::set<std::string> NameCont;

      /** \brief Create a data base object for the given FITS template file.
          \param tpl_file The name of the FITS template file used to create the file in memory.
      */
      PulsarDb(const std::string & tpl_file);

      /// \brief Destruct this PulsarDb object.
      virtual ~PulsarDb();

      /** \brief Load ephemerides and related information from the given ephemerides database file.
                 This copies the file in memory, and the version on disk will be unaffected.
          \param in_file The name of the input ephemerides database file.
      */
      virtual void load(const std::string & in_file);

      /** \brief Filter the current database using the given row filtering expression. The filtering is
                 performed in place.
          \param expression The row filtering expression. Examples: f2 != 0., #row > 50 && #row < 100
          \param filter_spin If true, spin parameters will be filtered by the given expression.
          \param filter_orbital If true, orbital parameters will be filtered by the given expression.
          \param filter_remark If true, ephemeris remarks will be filtered by the given expression.
      */
      virtual void filterExpression(const std::string & expression, bool filter_spin = true, bool filter_orbital = false,
        bool filter_remark = false);

      /** \brief Select ephemerides whose validation interval [VALID_SINCE, VALID_UNTIL] overlaps the given time range.
          \param t_start The start time of the interval.
          \param t_stop The stop time of the interval. (An exception will be thrown if t_stop < t_start).
      */
      virtual void filterInterval(double t_start, double t_stop);

      /** \brief Select ephemerides whose name matches the input name. Matching is case insensitive.

                 The name will be looked up two ways. It will be directly compared to the PSRNAME field in
                 the SPIN_PARAMETERS table. In addition, synonyms will be found in the ALTERNATIVE_NAMES table
                 and these will then be used to look up ephemerides from the SPIN_PARAMETERS table.
          \param pulsar_name The character string that contains the name of the pulsar. Examples: Crab, PSR B0531+21, PSR J0534+2200
      */
      virtual void filterName(const std::string & pulsar_name);

      /** \brief Select ephemerides whose SOLAR_SYSTEM_EPHEMERIS value matches the input solar_eph. Matching is case insensitive.
          \param solar_eph The name of the solar system ephemeris. Examples: JPL DE405, JPL DE200
      */
      virtual void filterSolarEph(const std::string & solar_eph);

      /** \brief Save the currently selected ephemerides into an output file.

                 All tables from the input file will be copied to the output file, but in each table, records which
                 do not match the current selection of ephemerides will be omitted. The matching is done based on PSRNAME for
                 the ORBITAL_PARAMETERS and ALTERNATIVE_NAMES, and on OBSERVER_CODE for the OBSERVERS table.
                 For example, if none of the ephemerides are for binary pulsars, the output ORBITAL_PARAMETERS
                 table will not include any ephemerides.
          \param out_file The name of the output file.
          \param creator The character string to be assigned to a value of CREATOR header keyword.
          \param author The character string to be assigned to a value of AUTHOR header keyword.
          \param clobber If true, it overwrites the output file even if it already exists.  If no and the output file
                 already exists, it throws an exception.
      */
      virtual void save(const std::string & out_file, const std::string & creator, const std::string & author,
        bool clobber = false) const;

      /** \brief Get the currently selected spin (pulsar) ephemerides.
          \param cont The container to fill the currently selected spin ephemerides in it. The previous contents of
                 the container will be removed by this method.
      */
      virtual void getEph(PulsarEphCont & cont) const { getEphBody(m_spin_par_table, m_spin_factory_cont, cont); }

      /** \brief Get the currently selected orbital ephemerides.
          \param cont The container to fill the currently selected orbital ephemerides in it. The spin ephemeris objects
                 that are stored in this container before calling this method will be destroyed and removed from the container.
      */
      virtual void getEph(OrbitalEphCont & cont) const { getEphBody(m_orbital_par_table, m_orbital_factory_cont, cont); }

      /** \brief Get the number of currently selected ephemerides.
          \param spin_table If true (default), the number of spin ephemerides is returned. If false, the number of orbital
                 ephemerides is returned.
      */
      virtual int getNumEph(bool spin_table = true) const;

      /** \brief Associate a keyword with a PulsarEph class as a handler of a FITS extension.
          \param eph_style Keyword to be associated to a PulsarEph class.
                 Note: Template parameter, EPHSTYLE, must be one of PulsarEph subclasses.
      */
      template <typename EPHSTYLE>
      void registerPulsarEph(const std::string & eph_style) {
        m_spin_factory_cont[eph_style] = &EphFactory<PulsarEph, EPHSTYLE>::getFactory();
      }

      /** \brief Associate a keyword with a OrbitalEph class as a handler of a FITS extension.
          \param eph_style Keyword to be associated to a OrbitalEph class.
                 Note: Template parameter, EPHSTYLE, must be one of OrbitalEph subclasses.
      */
      template <typename EPHSTYLE>
      void registerOrbitalEph(const std::string & eph_style) {
        m_orbital_factory_cont[eph_style] = &EphFactory<OrbitalEph, EPHSTYLE>::getFactory();
      }

      /** \brief Get the remarks stored in this pulsar ephemerides database.
          \param cont The container to fill the remarks in it. The EphStatus objects that are stored in this container
                 before calling this method will be destroyed and removed from the container.
      */
      virtual void getRemark(EphStatusCont & cont) const;

      /** \brief Get the creation history of this pulsar ephemerides database.
          \param command_history The container to fill the history of ephemeris loading and filtering.
          \param ancestry_record The container to fill the history records of all ephemerides database that have been loaded
                 to this ephemerides database.
      */
      virtual void getHistory(std::list<std::string> & command_history, std::list<std::string> & ancestry_record) const;

      /** \brief Examine the database to determine whether the named pulsar is in a binary system, and return the result.

                 This method determines that the database indicates the named pulsar is a binary pulsar if:
                 1) BINARY_FLAG column in SPIN_PARAMETERS extension for the named pulsar contains a logical true for
                    at least one entry for the named pulsar, or
                 2) an orbital ephemeris for the named pulsar is in ORBITAL_PARAMETERS extension.

                 The method determines that the database indicates the named pulsar is not a binary pulsar if:
                 3) BINARY_FLAG column in SPIN_PARAMETERS extension for the named pulsar contains a logical false for
                    at least one entry for the named pulsar.

                 The method returns true if either of conditions 1 or 2 is met AND condition 3 is not met. The method returns
                 false if none of conditions 1 and 2 is met AND condition 3 is met. The method returns false if none of the
                 above conditions is met. The method throws an exception if the above determinations conflict, i.e., if either
                 of conditions 1 and 2 is met at the same time as condition 3 is met.
          \param pulsar_name The character string that contains the name of the pulsar. Examples: Crab, PSR B0531+21, PSR J0534+2200
      */
      virtual bool isBinary(const std::string & pulsar_name) const;

      virtual void removeInterimFile();

    private:
      /** \brief Load ephemerides and related information from the given FITS file.
          \param in_file The name of the input FITS file.
      */
      virtual void loadFits(const std::string & in_file);

      /** \brief Helper method to find an extension that contains all required keyword-value pairs,
          and update its header if necessary.
          \param required_keyword List of keyword-value pairs that must be found in an extension to find.
      */
      virtual tip::Table * updateMatchingHeader(const tip::Header::KeySeq_t & required_keyword);

      /** \brief Load ephemerides and related information from the given TEXT file.
          \param in_file The name of the input TEXT file.
      */
      virtual void loadText(const std::string & in_file);

      /** \brief Clean up all extensions based on current set of selected spin and orbital ephemerides. All
          information in the OBSERVERS and ALTERNATIVE_NAMES extension which is not associated with a pulsar
          contained in the SPIN_PARAMETERS or ORBITAL_PARAMETERS extensions will be removed.
      */
      virtual void clean();

      /** \brief Creates a filtering expression from an input field name and a set of accepted values.
          \param field_name The name of the field on which to filter.
          \param values The set of allowed values.
      */
      virtual std::string createFilter(const std::string & field_name, const NameCont & values) const;

      /** \brief Helper method for getEph methods for PulsarEphCont and OrbitalEphCont.
          \param table_cont The container of tip tables, from which ephemerides are extracted.
          \param factory_cont The container of ephemeris factories to be used to create appropriate ephemeris objects.
          \param eph_cont The ephemeris container to be filled with appropriate ephemeris objects.
          Note: the original contents of this object will be removed from the container.
      */
      template <typename FactoryCont, typename EphCont>
      void getEphBody(const TableCont & table_cont, const FactoryCont & factory_cont, EphCont & eph_cont) const;

      /** \brief Helper method to collect all alternative names for the named pulsar.
          \param pulsar_name The character string that contains the name of the pulsar. Examples: Crab, PSR B0531+21, PSR J0534+2200
          \param alt_name_cont Container of alternative names for the named pulsar.
      */
      virtual void getAltName(const std::string & pulsar_name, NameCont & alt_name_cont) const;

      std::string m_tpl_file;
      tip::TipFile m_tip_file;
      TableCont m_all_table;
      TableCont m_spin_par_table;
      TableCont m_orbital_par_table;
      TableCont m_eph_remark_table;
      TableCont m_obs_code_table;
      TableCont m_psr_name_table;
      std::map<tip::Table *, int> m_table_generation_dict;
      std::map<std::string, IEphFactory<PulsarEph> *> m_spin_factory_cont;
      std::map<std::string, IEphFactory<OrbitalEph> *> m_orbital_factory_cont;
      std::list<std::string> m_command_history;
      std::list<std::string> m_ancestry_record;
  };

  template <typename FactoryCont, typename EphCont>
  void PulsarDb::getEphBody(const TableCont & table_cont, const FactoryCont & factory_cont, EphCont & eph_cont) const {
    // Empty container then refill it.
    for (typename EphCont::reverse_iterator itor = eph_cont.rbegin(); itor != eph_cont.rend(); ++itor) delete *itor;
    eph_cont.clear();

    // Reserve space for all ephemerides.
    int num_record = 0;
    for (TableCont::const_iterator itor = table_cont.begin(); itor != table_cont.end(); ++itor) {
      num_record += (*itor)->getNumRecords();
    }
    eph_cont.reserve(num_record);

    for (TableCont::const_iterator table_itor = table_cont.begin(); table_itor != table_cont.end(); ++table_itor) {
      const tip::Table & table = **table_itor;
      const tip::Header & header(table.getHeader());

      // Get the extension name.
      std::string ext_name;
      header["EXTNAME"].get(ext_name);

      // Check and read EPHSTYLE keyword to select a proper ephemeris factory.
      if (header.find("EPHSTYLE") == header.end()) {
        // Note: EPHSTYLE must exist in SPIN_PARAMETERS and ORBITAL_PARAMETERS extensions, and it is enforced
        //       in the constructor of this class. Not finding EPHSTYLE here suggests inconsistency in methods
        //       of this class, or a bug most likely.
        throw std::logic_error("EPHSTYLE header keyword is missing in " + ext_name + " extension");
      }
      std::string eph_style;
      header["EPHSTYLE"].get(eph_style);

      // Use a registered subclass of PulsarEph or OrbitalEph whichever appropriate, if EPHSTYLE keyword exists.
      typename FactoryCont::mapped_type factory(0);
      typename FactoryCont::const_iterator factory_itor = factory_cont.find(eph_style);
      if (factory_itor != factory_cont.end()) {
        factory = factory_itor->second;
      } else {
        throw std::runtime_error("Unknown ephemeris style for " + ext_name + " extension: EPHSTYLE = " + eph_style);
      }

      // Iterate over current selection.
      for (tip::Table::ConstIterator record_itor = table.begin(); record_itor != table.end(); ++record_itor) {
        // For convenience, get record from iterator.
        tip::Table::ConstRecord & record(*record_itor);

        // Add the ephemeris to the container.
        eph_cont.push_back(factory->create(record, header));
      }
    }
  }

}

#endif
