/** \file PeriodicityTestApp.cxx
    \brief Implmentation of PeriodicityTestApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PeriodicityTestApp.h"

#include <ctime>
#include <limits>
#include <list>
#include <memory>
#include <stdexcept>

#include "facilities/commonUtilities.h"

#include "hoops/hoops.h"

#include "st_app/AppParGroup.h"

#include "st_facilities/FileSys.h"

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "ChiSquaredTestArray.h"
#include "HTestArray.h"
#include "PeriodicityTestArray.h"
#include "RayleighTestArray.h"
#include "Z2nTestArray.h"

static const std::string s_cvs_id = "$Name: ScienceTools-v10r0p5-fssc-20150518 $";

PeriodicityTestApp::PeriodicityTestApp(): m_os("PeriodicityTestApp", "", 2) {
  setName("gtptest");
  setVersion(s_cvs_id);
}

PeriodicityTestApp::~PeriodicityTestApp() throw() {}

void PeriodicityTestApp::run() {
  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup());

  // Prompt and save.
  pars.Prompt();
  pars.Save();

  // Read parameters to open input event files.
  std::string event_file = pars["evfile"];
  std::string event_extension = pars["evtable"];
  std::string field_name = pars["pphasefield"];

  // Open the event table(s) for reading.
  typedef std::list<const tip::Table *> table_list_type;
  table_list_type table_list;
  st_facilities::FileSys::FileNameCont file_name_cont = st_facilities::FileSys::expandFileList(event_file);
  for (st_facilities::FileSys::FileNameCont::const_iterator itor = file_name_cont.begin(); itor != file_name_cont.end(); ++itor) {
    std::string file_name = *itor;
    const tip::Table * table = tip::IFileSvc::instance().readTable(file_name, event_extension);

    // Check whether pulse phase field exists.
    try {
      table->getFieldIndex(field_name);
    } catch (const tip::TipException &) {
      throw std::runtime_error("Input file \"" + file_name + "\" does not contain a pulse phase column named \"" + field_name + "\"");
    }

    // Add it to the table list.
    table_list.push_back(table);
  }

  // Create a list of periodicity tests.
  typedef std::list<PeriodicityTestArray *> test_list_type;
  test_list_type test_list;
  std::map<std::string, std::string> ext_name_dict;

  // Add periodicity tests to the list.
  long num_phase = pars["numphase"];
  PeriodicityTestArray * chi2_test_array(new ChiSquaredTestArray(1, num_phase));
  test_list.push_back(chi2_test_array);
  ext_name_dict[chi2_test_array->getTestName()] = "CHI2TEST";

  PeriodicityTestArray * rayleigh_test_array(new RayleighTestArray(1));
  test_list.push_back(rayleigh_test_array);
  ext_name_dict[rayleigh_test_array->getTestName()] = "RAYLEIGHTEST";

  long num_harm = pars["numharm"];
  PeriodicityTestArray * z2n_test_array(new Z2nTestArray(1, num_harm));
  test_list.push_back(z2n_test_array);
  ext_name_dict[z2n_test_array->getTestName()] = "Z2NTEST";

  long max_harm = pars["maxharm"];
  PeriodicityTestArray * h_test_array(new HTestArray(1, max_harm));
  test_list.push_back(h_test_array);
  ext_name_dict[h_test_array->getTestName()] = "HTEST";

  // Loop over events in the event table(s).
  for (table_list_type::const_iterator table_itor = table_list.begin(); table_itor != table_list.end(); ++table_itor) {
    const tip::Table & table = **table_itor;
    for (tip::Table::ConstIterator event_itor = table.begin(); event_itor != table.end(); ++event_itor) {
      const tip::ConstTableRecord & record = *event_itor;

      // Read phase value from this event, as a signle variable of double type.
      double phase_value;
      record[field_name].get(phase_value);

      // Fill the phase value into the tests.
      for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) {
        PeriodicityTestArray & test_array = **itor;
        test_array.fill(0, phase_value);
      }
    }
  }

  // Update the data contents of viewer.
  for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) {
    PeriodicityTestArray & test_array = **itor;
    test_array.updateViewer(0);
  }

  // Get parameters for output.
  std::string out_file = pars["outfile"];
  bool plot = pars["plot"];
  std::string title = pars["title"];
  bool clobber = pars["clobber"];

  // Use default title if user did not specify one.
  std::string title_uc(title);
  for (std::string::iterator itor = title_uc.begin(); itor != title_uc.end(); ++itor) *itor = std::toupper(*itor);
  if (title_uc != "DEFAULT") {
    for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) {
      PeriodicityTestArray & test_array = **itor;
      test_array.getViewer().setTitle(title);
    }
  }

  // Interpret output file parameter.
  std::string out_file_uc = out_file;
  for (std::string::iterator itor = out_file_uc.begin(); itor != out_file_uc.end(); ++itor) *itor = std::toupper(*itor);

  if ("NONE" != out_file_uc) {
    // Find the template file.
    using namespace facilities;
    std::string template_file = commonUtilities::joinPath(commonUtilities::getDataPath("periodSearch"), "periodicity-test-out.tpl");

    // Create output file.
    tip::IFileSvc::instance().createFile(out_file, template_file, clobber);

    // Construct a character string representing file creation time in UTC.
    // Note: UTC is the default time system for DATE header keyword in the FITS standard.
    std::time_t current_time = std::time(0);
    struct std::tm * gm_time_struct = std::gmtime(&current_time);
    char gm_time_char[] = "YYYY-MM-DDThh:mm:ss";
    std::strftime(gm_time_char, sizeof(gm_time_char), "%Y-%m-%dT%H:%M:%S", gm_time_struct);

    // Create a header line for HISTORY records.
    std::string creator_name = getName() + " " + getVersion();
    std::string file_creation_time(gm_time_char);
    std::string header_line("File created by " + creator_name + " on " + file_creation_time);

    // Loop over periodicity tests.
    for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) {
      PeriodicityTestArray & test_array = **itor;

      // Open output file.
      std::string ext_name = ext_name_dict[test_array.getTestName()];
      std::auto_ptr<tip::Table> out_table(tip::IFileSvc::instance().editTable(out_file, ext_name));

      // Write the summary to the output header, and the data to the output table.
      test_array.getViewer().write(*out_table);

      // Update header keywords.
      tip::Header & header(out_table->getHeader());
      tip::Header::KeyValCont_t keywords;
      keywords.push_back(tip::Header::KeyValPair_t("DATE", file_creation_time));
      keywords.push_back(tip::Header::KeyValPair_t("CREATOR", creator_name));
      keywords.push_back(tip::Header::KeyValPair_t("DATASUM", "-1")); // Force update of DATASUM keyword.
      header.update(keywords);

      // Write out all the parameters into HISTORY keywords.
      header.addHistory(header_line);
      const st_app::AppParGroup & const_pars(pars);
      for (hoops::ConstGenParItor par_itor = const_pars.begin(); par_itor != const_pars.end(); ++par_itor) {
        std::ostringstream oss_par;
        oss_par << getName() << ".par: " << **par_itor;
        header.addHistory(oss_par.str());
      }
    }
  }

  // Compute the statistical test results and write them to the screen.
  for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) {
    PeriodicityTestArray & test_array = **itor;
    test_array.getViewer().write(m_os, 2, 5);
  }

  // Display a plot, if desired.
  if (plot) {
    StatisticViewer & viewer = chi2_test_array->getViewer();
    viewer.setLabel(0, "Pulse Phase");
    viewer.setLabel(1, "Counts");
    viewer.plot();
  }

  // Delete the event table(s).
  for (table_list_type::reverse_iterator itor = table_list.rbegin(); itor != table_list.rend(); ++itor) delete *itor;

  // Delete the peridicity tests.
  for (test_list_type::reverse_iterator itor = test_list.rbegin(); itor != test_list.rend(); ++itor) delete *itor;
}
