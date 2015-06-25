/** \file BinConfig.cxx
    \brief Implementation of helper class which uses standard sets of parameters to configure binners for standard applications.
    \author James Peachey, HEASARC
*/
#include <algorithm>
#include <cctype>
#include <climits>
#include <memory>
#include <stdexcept>
#include <string>
#include <cmath>

#include "GlastGbmBinConfig.h"
#include "GlastLatBinConfig.h"
#include "evtbin/BinConfig.h"
#include "evtbin/ConstSnBinner.h"
#include "evtbin/Gti.h"
#include "evtbin/LinearBinner.h"
#include "evtbin/LogBinner.h"
#include "evtbin/OrderedBinner.h"

// Interactive parameter file access from st_app.
#include "st_app/AppParGroup.h"

#include "st_facilities/FileSys.h"

#include "tip/Extension.h"
#include "tip/FileSummary.h"
#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

namespace evtbin {

  BinConfig::ConfigCont BinConfig::s_config_cont;

  void BinConfig::load() {
    GlastGbmBinConfig::load();
    GlastLatBinConfig::load();
  }

  BinConfig * BinConfig::create(const std::string & ev_file_name) {
    using namespace st_facilities;
    FileSys::FileNameCont file_name_cont = FileSys::expandFileList(ev_file_name);
    BinConfig * config = 0;
    std::string mission;
    std::string instrument;
    for (FileSys::FileNameCont::iterator file_itor = file_name_cont.begin(); file_itor != file_name_cont.end(); ++file_itor) {
      // Get container of all extensions in file.
      tip::FileSummary summary;
      tip::IFileSvc::instance().getFileSummary(*file_itor, summary);
  
      // Read all extensions in input file, looking for mission and instrument names.
      for (tip::FileSummary::const_iterator itor = summary.begin(); itor != summary.end(); ++itor) {
        // Open extension.
        std::auto_ptr<const tip::Extension> table(tip::IFileSvc::instance().readExtension(*file_itor, itor->getExtId()));
  
        try {
          if (mission.empty()) table->getHeader()["TELESCOP"].get(mission);
          if (instrument.empty()) table->getHeader()["INSTRUME"].get(instrument);
          break;
        } catch (const tip::TipException &) {
          // Read each extension until these keywords are found.
          continue;
        }
      }
    }

    if (!mission.empty() && !instrument.empty()) {
      // Put mission name in all uppercase to deal with any mixed cases
      std::transform(mission.begin(),mission.end(),mission.begin(),::toupper);
      // Find a prototype for a bin configuration appropriate for this mission/instrument.
      ConfigCont::iterator found = s_config_cont.find(mission + "::" + instrument);
      if (s_config_cont.end() != found) {
        config = found->second->clone();
      } else {
        throw std::runtime_error("BinConfig::create was unable to find a configuration for mission \"" + mission +
          "\", instrument \"" + instrument + "\" while processing file \"" +
          ev_file_name + "\"");
      }
    }

    if (0 == config)
      throw std::runtime_error("BinConfig::create was unable to determine the mission/instrument in file \"" +
        ev_file_name + "\"");

    return config;
  }

  BinConfig::~BinConfig() {
    // This object is going away. Remove it from the container of prototypes if it is there.
    for (ConfigCont::iterator itor = s_config_cont.begin(); itor != s_config_cont.end(); ++itor) {
      if (this == itor->second) {
        s_config_cont.erase(itor);
        break;
      }
    }
  }

  void BinConfig::energyParPrompt(st_app::AppParGroup & par_group) const {
    parPrompt(par_group, "ebinalg", "efield", "emin", "emax", "denergy", "enumbins", "ebinfile");
  }

  void BinConfig::spatialParPrompt(st_app::AppParGroup & par_group) const {
    par_group.Prompt("nxpix");
    par_group.Prompt("nypix");
    par_group.Prompt("binsz");
    par_group.Prompt("coordsys");
    par_group.Prompt("xref");
    par_group.Prompt("yref");
    par_group.Prompt("rafield");
    par_group.Prompt("decfield");
    par_group.Prompt("axisrot");
    par_group.Prompt("proj");
  }

  void BinConfig::healpixParPrompt(st_app::AppParGroup & par_group) const {
    par_group.Prompt("hpx_ordering_scheme");
    par_group.Prompt("hpx_order");
    par_group.Prompt("coordsys");
    par_group.Prompt("hpx_ebin");
   
  }

  void BinConfig::timeParPrompt(st_app::AppParGroup & par_group) const {
    parPrompt(par_group, "tbinalg", "tfield", "tstart", "tstop", "dtime", "ntimebins", "tbinfile", "snratio",
      "lcemin", "lcemax");
  }

  Binner * BinConfig::createEnergyBinner(const st_app::AppParGroup & par_group) const {
    return createBinner(par_group, "ebinalg", "efield", "emin", "emax", "denergy", "enumbins", "ebinfile",
      "ENERGYBINS", "E_MIN", "E_MAX");
  }

  Binner * BinConfig::createTimeBinner(const st_app::AppParGroup & par_group) const {
    return createBinner(par_group, "tbinalg", "tfield", "tstart", "tstop", "dtime", "ntimebins", "tbinfile",
      "TIMEBINS", "START", "STOP", "snratio", "lcemin", "lcemax");
  }

  Binner * BinConfig::createEbounds(const st_app::AppParGroup & par_group) const {
    return createBinner(par_group, "ebinalg", "efield", "emin", "emax", "denergy", "enumbins", "ebinfile",
      "ENERGYBINS", "E_MIN", "E_MAX");
  }

  Gti * BinConfig::createGti(const st_app::AppParGroup & par_group) const {
    return new Gti(par_group["evfile"]);
  }

  void BinConfig::timeParDefaults(st_app::AppParGroup & par_group, const std::string & timevalue) const {
    // In the case of INDEF, UNDEF, or a null value, get value from the original fits header.
    double headerTime;
    std::string timevalUpper;
    timevalUpper=timevalue;
    for (std::string::iterator itor = timevalUpper.begin(); itor != timevalUpper.end(); ++itor) *itor = toupper(*itor);
    try {
      par_group.Prompt(timevalue);
    }catch(const hoops::Hexception &){
      std::auto_ptr<const tip::Extension> table(tip::IFileSvc::instance().readExtension(par_group["evfile"].Value(), "EVENTS"));
      table->getHeader()[timevalUpper].get(headerTime);
      par_group[timevalue]=headerTime;
    }
  }

  void BinConfig::parPrompt(st_app::AppParGroup & par_group, const std::string & alg, const std::string & in_field,
    const std::string & bin_begin, const std::string & bin_end, const std::string & bin_size, const std::string & num_bins,
    const std::string & bin_file, const std::string & sn_ratio, const std::string & lc_emin, const std::string & lc_emax) const {
    // Determine the time binning algorithm.
    par_group.Prompt(alg);
    par_group.Prompt(in_field);

    // Get the type of bin specification.
    std::string bin_type = par_group[alg].Value();

    // Make all upper case for case-insensitive comparisons.
    for (std::string::iterator itor = bin_type.begin(); itor != bin_type.end(); ++itor) *itor = toupper(*itor);

    if (bin_type == "LIN") {
      // Get remaining parameters needed for linearly uniform interval binner.
      if (bin_begin == "tstart") {
	timeParDefaults(par_group,bin_begin);
	timeParDefaults(par_group,bin_end);
      } else {
	par_group.Prompt(bin_begin);
	par_group.Prompt(bin_end);
      }
      par_group.Prompt(bin_size);
      double begin = par_group[bin_begin];
      double end = par_group[bin_end];
      double size = par_group[bin_size];
      if(size == 0) throw std::runtime_error("Bins must have non-zero width!");
      if(ceil((end - begin)/size) > LONG_MAX) throw std::length_error("Number of bins exceeds max possible on this system!");
    } else if (bin_type == "LOG") {
      // Get remaining parameters needed for logarithmically uniform interval binner.
      par_group.Prompt(bin_begin);
      par_group.Prompt(bin_end);
      par_group.Prompt(num_bins);
    } else if (bin_type == "FILE") {
      // Get remaining parameters needed for user defined bins from a bin file.
      par_group.Prompt(bin_file);
    } else if (bin_type == "SNR") {
      if (bin_begin == "tstart") {
	timeParDefaults(par_group,bin_begin);
	timeParDefaults(par_group,bin_end);
      } else {
	par_group.Prompt(bin_begin);
	par_group.Prompt(bin_end);
      }
      par_group.Prompt(sn_ratio);
      par_group.Prompt(lc_emin);
      par_group.Prompt(lc_emax);
    } else throw std::runtime_error(std::string("Unknown binning algorithm ") + par_group[alg].Value());
  }

  Binner * BinConfig::createBinner(const st_app::AppParGroup & par_group, const std::string & alg,
    const std::string & in_field, const std::string & bin_begin, const std::string & bin_end, const std::string & bin_size,
    const std::string & num_bins, const std::string & bin_file, const std::string & bin_ext, const std::string & start_field,
    const std::string & stop_field, const std::string & sn_ratio, const std::string & lc_emin, const std::string & lc_emax) const {
    using namespace evtbin;

    Binner * binner = 0;

    // Get the type of bin specification.
    std::string bin_type = par_group[alg].Value();

    // Make all upper case for case-insensitive comparisons.
    for (std::string::iterator itor = bin_type.begin(); itor != bin_type.end(); ++itor) *itor = toupper(*itor);

    if (bin_type == "LIN") {
      binner = new LinearBinner(par_group[bin_begin], par_group[bin_end], par_group[bin_size], par_group[in_field]);
    } else if (bin_type == "LOG") {
      binner = new LogBinner(par_group[bin_begin], par_group[bin_end], par_group[num_bins], par_group[in_field]);
    } else if (bin_type == "FILE") {
      // Create interval container for user defined bin intervals.
      OrderedBinner::IntervalCont_t intervals;

      std::string bin_ext_uc = bin_ext;

      // Make all upper case for case-insensitive comparisons.
      for (std::string::iterator itor = bin_ext_uc.begin(); itor != bin_ext_uc.end(); ++itor) *itor = toupper(*itor);

      // Open the data file.
      std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(par_group[bin_file], bin_ext));

// TODO Refactor this!
      // Temporary hack pending a complete handling of units in tip.
      double factor = 1.;

      // For the moment, energy bins are in keV, period, so in the case of an energy binner, convert keV to MeV.
      if (bin_ext_uc == "EBOUNDS" || bin_ext_uc == "ENERGYBINS") factor = 1.e-3;

      // Iterate over the file, saving the relevant values into the interval array.
      for (tip::Table::ConstIterator itor = table->begin(); itor != table->end(); ++itor) {
        intervals.push_back(Binner::Interval(factor * (*itor)[start_field].get(), factor * (*itor)[stop_field].get()));
      }

      // Create binner from these intervals.
      binner = new OrderedBinner(intervals, par_group[in_field]);

    } else if (bin_type == "SNR") {
      binner = new ConstSnBinner(par_group[bin_begin], par_group[bin_end], par_group[sn_ratio], par_group[lc_emin],
        par_group[lc_emax], std::vector<double>(), par_group[in_field]);
    } else throw std::runtime_error(std::string("Unknown binning algorithm ") + par_group[alg].Value());

    return binner;
  }

}
