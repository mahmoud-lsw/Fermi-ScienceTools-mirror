/** \file PowerSpectrumApp.cxx
    \brief Implmentation of PowerSpectrumApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PowerSpectrumApp.h"

#include <cctype>
#include <ctime>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/commonUtilities.h"

#include "hoops/hoops.h"

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphStatus.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeSystem.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "FourierAnalysis.h"
#include "PeriodSearch.h"
#include "StatisticViewer.h"

static const std::string s_cvs_id = "$Name: ScienceTools-v10r0p5-fssc-20150518 $";

PowerSpectrumApp::PowerSpectrumApp(): m_os("PowerSpectrumApp", "", 2) {
  setName("gtpspec");
  setVersion(s_cvs_id);

  st_app::AppParGroup & pars(getParGroup());
  pars.setSwitch("ephstyle");
  pars.setCase("ephstyle", "FREQ", "f1f0ratio");
  pars.setCase("ephstyle", "FREQ", "f2f0ratio");
  pars.setCase("ephstyle", "PER", "p1p0ratio");
  pars.setCase("ephstyle", "PER", "p2p0ratio");
}

PowerSpectrumApp::~PowerSpectrumApp() throw() {}

void PowerSpectrumApp::runApp() {
  m_os.setMethod("runApp()");
  st_app::AppParGroup & pars(getParGroup());

  // Prompt for all parameters and save them.
  prompt(pars);

  // Get parameters.
  std::string out_file = pars["outfile"];
  long num_bins = pars["numbins"];
  double low_f_cut = pars["lowfcut"];
  bool plot = pars["plot"];
  std::string title = pars["title"];
  bool clobber = pars["clobber"];

  // Open the event file(s).
  openEventFile(pars);

  // Handle leap seconds.
  std::string leap_sec_file = pars["leapsecfile"];
  timeSystem::TimeSystem::setDefaultLeapSecFileName(leap_sec_file);

  // Setup time correction mode.
  defineTimeCorrectionMode("NONE", SUPPRESSED, SUPPRESSED, SUPPRESSED);
  defineTimeCorrectionMode("AUTO", ALLOWED,    ALLOWED,    ALLOWED);
  defineTimeCorrectionMode("BARY", REQUIRED,   SUPPRESSED, SUPPRESSED);
  defineTimeCorrectionMode("BIN",  REQUIRED,   REQUIRED,   SUPPRESSED);
  defineTimeCorrectionMode("PDOT", REQUIRED,   SUPPRESSED, REQUIRED);
  defineTimeCorrectionMode("ALL",  REQUIRED,   REQUIRED,   REQUIRED);
  selectTimeCorrectionMode(pars);

  // Set up EphComputer for arrival time corrections.
  pulsarDb::SloppyEphChooser chooser;
  initEphComputer(pars, chooser, "NONE", m_os.info(4));

  // Use user input (parameters) together with computer to determine corrections to apply.
  bool vary_ra_dec = false;
  bool guess_pdot = false;
  initTimeCorrection(pars, vary_ra_dec, guess_pdot, m_os.info(3));

  // Report ephemeris status.
  std::set<pulsarDb::EphStatusCodeType> code_to_report;
  code_to_report.insert(pulsarDb::Remarked);
  reportEphStatus(m_os.warn(), code_to_report);

  // Compute start time of the data set.
  double tstart = computeElapsedSecond(getStartTime());
  double tstop = computeElapsedSecond(getStopTime());

  // Get binwidth parameter.
  double bin_width = pars["binwidth"];

  // Create the proper period search object..
  std::auto_ptr<PeriodSearch> search(new FourierAnalysis(tstart, tstop, bin_width, num_bins, "Hz"));

  for (setFirstEvent(); !isEndOfEventList(); setNextEvent()) {
    // Get event time as AbsoluteTime.
    timeSystem::AbsoluteTime abs_evt_time(getEventTime());

    // Convert event time to target time representation.
    double target_evt_time = computeElapsedSecond(abs_evt_time);

    // Fill into the period search object.
    search->fill(target_evt_time);
  }

  // Compute the statistics.
  search->computeStat();

  // Update the statistic viewer in PeriodSearch object, and get a reference to it.
  search->updateViewer(low_f_cut);
  StatisticViewer & viewer(search->getViewer());

  // Set a plot title: use default title if user did not specify one.
  std::string title_uc(title);
  for (std::string::iterator itor = title_uc.begin(); itor != title_uc.end(); ++itor) *itor = std::toupper(*itor);
  if (title_uc != "DEFAULT") viewer.setTitle(title);

  // Set unit for a plot.
  viewer.setUnit(0, "Hz");

  // Interpret output file parameter.
  std::string out_file_uc = out_file;
  for (std::string::iterator itor = out_file_uc.begin(); itor != out_file_uc.end(); ++itor) *itor = std::toupper(*itor);

  if ("NONE" != out_file_uc) {
    // Find the template file.
    using namespace facilities;
    std::string template_file = commonUtilities::joinPath(commonUtilities::getDataPath("periodSearch"), "period-search-out.tpl");

    // Create output file.
    tip::IFileSvc::instance().createFile(out_file, template_file, clobber);

    // Create a header line for HISTORY records.
    std::string creator_name = getName() + " " + getVersion();
    std::string file_creation_time(createUtcTimeString());
    std::string header_line("File created by " + creator_name + " on " + file_creation_time);

    // Open output file.
    std::auto_ptr<tip::Table> out_table(tip::IFileSvc::instance().editTable(out_file, "POWER_SPECTRUM"));

    // Write the summary to the output header, and the data to the output table.
    viewer.write(*out_table);

    // Get a target name to put in the output file.
    std::string psr_name = pars["psrname"];
    for (std::string::iterator itor = psr_name.begin(); itor != psr_name.end(); ++itor) *itor = std::toupper(*itor);

    // Update header keywords.
    tip::Header & header(out_table->getHeader());
    tip::Header::KeyValCont_t keywords;
    keywords.push_back(tip::Header::KeyValPair_t("DATE", file_creation_time));
    keywords.push_back(tip::Header::KeyValPair_t("CREATOR", creator_name));
    keywords.push_back(tip::Header::KeyValPair_t("OBJECT", psr_name));
    keywords.push_back(tip::Header::KeyValPair_t("DATASUM", "-1")); // Force update of DATASUM keyword.
    header.update(keywords);

    // Write out all the parameters into HISTORY keywords.
    writeParameter(pars, header_line, header);
  }

  // Write the search results to the screen.
  viewer.write(m_os, 2, 5);

  // Display a plot, if desired.
  if (plot) {
    viewer.setLabel(0, "Pulse Frequency");
    viewer.setLabel(1, "Fourier Power");
    viewer.plot();
  }
}

void PowerSpectrumApp::prompt(st_app::AppParGroup & pars) {
  // Prompt for most parameters automatically.
  pars.Prompt("evfile");
  pars.Prompt("evtable");
  pars.Prompt("timefield");
  pars.Prompt("scfile");
  pars.Prompt("sctable");
  pars.Prompt("psrdbfile");
  pars.Prompt("psrname");
  pars.Prompt("outfile");
  pars.Prompt("binwidth");
  pars.Prompt("numbins");
  pars.Prompt("lowfcut");

  pars.Prompt("timeorigin");
  std::string origin_style = pars["timeorigin"];
  for (std::string::iterator itor = origin_style.begin(); itor != origin_style.end(); ++itor) *itor = std::toupper(*itor);
  if (origin_style == "USER") {
    pars.Prompt("usertime");
    pars.Prompt("userformat");
    pars.Prompt("usersys");
  }

  pars.Prompt("ra");
  pars.Prompt("dec");

  // Prompt for f1 & f2 / p1 & p2 even if pdot correction is NOT selected.
  pars.Prompt("ephstyle");
  std::string eph_style = pars["ephstyle"];
  for (std::string::iterator itor = eph_style.begin(); itor != eph_style.end(); ++itor) *itor = std::toupper(*itor);
  if (eph_style == "FREQ") {
    pars.Prompt("f1f0ratio");
    pars.Prompt("f2f0ratio");
  } else if (eph_style == "PER") {
    pars.Prompt("p1p0ratio");
    pars.Prompt("p2p0ratio");
  } else {
    throw std::runtime_error("Unknown ephemeris style " + eph_style);
  }

  pars.Prompt("tcorrect");
  pars.Prompt("solareph");
  pars.Prompt("matchsolareph");
  pars.Prompt("angtol");
  pars.Prompt("plot");
  pars.Prompt("title");
  pars.Prompt("leapsecfile");
  pars.Prompt("reportephstatus");
  pars.Prompt("chatter");
  pars.Prompt("clobber");
  pars.Prompt("debug");
  pars.Prompt("gui");
  pars.Prompt("mode");

  // Save current values of the parameters.
  pars.Save();
}
