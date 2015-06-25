/** \file test_evtbin.cxx
    \brief Event bin test program.
    \author Yasushi Ikebe, GSSC
            James Peachey, HEASARC
*/

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

// Class encapsulating a Bayesian block binner.
#include "evtbin/BayesianBinner.h"
// Class encapsulating a binner configuration helper object.
#include "evtbin/BinConfig.h"
// Class encapsulating a count map.
#include "evtbin/CountMap.h"
// Class encapsulating a count cube.
#include "evtbin/CountCube.h"
// Glass encapsulating GTIs
#include "evtbin/Gti.h"
// Class encapsulating a 1 dimensional histogram.
#include "evtbin/Hist1D.h"
// Class encapsulating a 2 dimensional histogram.
#include "evtbin/Hist2D.h"
// Light curve abstractions.
#include "evtbin/LightCurve.h"
// Class encapsulating description of a const s/n binner.
#include "evtbin/ConstSnBinner.h"
// Class encapsulating description of a binner with linear equal size bins.
#include "evtbin/LinearBinner.h"
// Class encapsulating description of a binner with logarithmically equal size bins.
#include "evtbin/LogBinner.h"
// Class encapsulating description of a binner with ordered but otherwise arbitrary bins.
#include "evtbin/OrderedBinner.h"
// Class encapsulating description of a HEALPIX binner 
#include "evtbin/HealpixBinner.h"
// Class encapsulating description of a HEALPIX map
#include "evtbin/HealpixMap.h"
// Class for binning into Hist objects from tip objects:
#include "evtbin/RecordBinFiller.h"
// Multiple spectra abstractions.
#include "evtbin/MultiSpec.h"
// Single spectrum abstractions.
#include "evtbin/SingleSpec.h"
// Application parameter class.
#include "st_app/AppParGroup.h"
// Application base class.
#include "st_app/StApp.h"
// Factory used by st_app's standard main to create application object.
#include "st_app/StAppFactory.h"
// Utility used to find test data for this application.
#include "st_facilities/Env.h"
// Utility used to find files for this application.
#include "st_facilities/FileSys.h"
// Message utilities.
#include "st_stream/st_stream.h"
#include "st_stream/StreamFormatter.h"
// Tip File access.
#include "tip/IFileSvc.h"
// Tip Table access.
#include "tip/Table.h"
// Tip type definitions
#include "tip/tip_types.h"

#include "facilities/commonUtilities.h"
using namespace evtbin;

const std::string s_cvs_id("$Name: ScienceTools-v10r0p5-fssc-20150518 $");

/** \class EvtBinTest
    \brief Application singleton for evtbin test program.
*/
class EvtBinTest : public st_app::StApp {
  public:
    EvtBinTest();

    virtual ~EvtBinTest() throw() {}

    /** \brief Run all tests.
    */
    virtual void run();

    void testLinearBinner();

    void testLogBinner();

    void testOrderedBinner();

    void testHealpixBinner();

    void testHealpixMap();

    void testHist1D();

    void testHist2D();

    void testLightCurve();

    void testSingleSpectrum();

    void testMultiSpectra();

    void testCountMap();

    void testCountCube();

    void testBinConfig();

    void testGti();

    void testConstSnBinner();

    void testBayesianBinner();

    void testMultipleFiles();

  private:
    st_stream::StreamFormatter m_os;
    std::string m_data_dir;
    std::string m_ft1_file;
    std::string m_ft2_file;
    std::string m_gbm_file;
    double m_t_start;
    double m_t_stop;
    double m_e_min;
    double m_e_max;
    double m_gbm_t_start;
    double m_gbm_t_stop;
    bool m_failed;
};

EvtBinTest::EvtBinTest(): m_os("EvtBinTest", "EvtBinTest", 2), m_t_start(223850000.0), m_t_stop(223950000.0),
  m_e_min(30.), m_e_max(6000.), m_gbm_t_start(-4.193924833089113E-03), m_gbm_t_stop(1.368369758129120E-01), m_failed(false) {
  setName("test_evtbin");
  setVersion(s_cvs_id);
  m_data_dir = facilities::commonUtilities::getDataPath("evtbin");
  m_ft1_file = facilities::commonUtilities::joinPath(m_data_dir, "ft1tiny.fits");
  m_ft2_file = facilities::commonUtilities::joinPath(m_data_dir, "ft2tiny.fits");
  m_gbm_file = facilities::commonUtilities::joinPath(m_data_dir, "gbmtiny.fits");
}

void EvtBinTest::run() {
  m_failed = false;

  std::cerr.precision(std::numeric_limits<double>::digits10);
  std::cout.precision(std::numeric_limits<double>::digits10);

  // Test high level bin configuration object first, as it tests some things which need to be done first.
  testBinConfig();
  // Test linear binner:
  testLinearBinner();
  // Test logarithmic binner:
  testLogBinner();
  // Test ordered binner:
  testOrderedBinner();
  // Test healpix binner:
  testHealpixBinner();
  // Test healpix map:
  testHealpixMap();
  // Test one dimensional histogram:
  testHist1D();
  // Test two dimensional histogram:
  testHist2D();
  // Test light curve with no energy binning (using Tip):
  testLightCurve();
  // Test single spectrum with no time binning (using Tip):
  testSingleSpectrum();
  // Test multiple spectra with time binning (using Tip):
  testMultiSpectra();
  // Test count map (using Tip):
  testCountMap();
  // Test count cube (using Tip):
  testCountCube();
  // Test high level bin configuration object.
  testGti();
  // Test const s/n binner:
  testConstSnBinner();
  // Test Bayesian block binner:
  testBayesianBinner();
  // Test getting input from multiple files:
  testMultipleFiles();

  // Report problems, if any.
  if (m_failed) throw std::runtime_error("Unit test failed");
}

void EvtBinTest::testLinearBinner() {
  std::string msg;

  // Create a linear binner with bin width == 15. spanning the interval [0, 100):
  LinearBinner binner(0., 100., 15.);

  // Make sure there are 7 bins:
  msg = "LinearBinner::getNumBins()";
  if (7 != binner.getNumBins()) {
    std::cerr << msg << " returned " << binner.getNumBins() << ", not 7" << std::endl;
    m_failed = true;
  }

  // Make sure values are correctly classified:
  msg = "LinearBinner::computeIndex(";
  for (long ii = 0; ii < 7; ++ii) {
    int value = 9 + ii * 15;
    long index = binner.computeIndex(value);
    if (index != ii) {
      std::cerr << msg << value << ") returned " << index << ", not " << ii << std::endl;
      m_failed = true;
    }
  }

  // Left endpoint should be included:
  long index = binner.computeIndex(0);
  if (0 != index) {
    std::cerr << msg << "0) returned " << index << ", not 0" << std::endl;
    m_failed = true;
  }

  // Right endpoint should be excluded:
  index = binner.computeIndex(100);
  if (-1 != index) {
    std::cerr << msg << "100) returned " << index << ", not -1" << std::endl;
    m_failed = true;
  }

  // Left of left endpoint should return index < 0:
  index = binner.computeIndex(-1);
  if (0 <= index) {
    std::cerr << msg << "-1) returned " << index << ", which is >= 0" << std::endl;
    m_failed = true;
  }

  // Right of right endpoint should return index < 0:
  index = binner.computeIndex(101);
  if (0 <= index) {
    std::cerr << msg << "101) returned " << index << ", which is >= 0" << std::endl;
    m_failed = true;
  }

  // Make sure nice symmetric intervals are also handled correctly.
  LinearBinner binner2(0., 100., 10.);

  // Make sure there are 10 bins:
  msg = "LinearBinner::getNumBins()";
  if (10 != binner2.getNumBins()) {
    std::cerr << msg << " returned " << binner.getNumBins() << ", not 10" << std::endl;
    m_failed = true;
  }

}

void EvtBinTest::testLogBinner() {
  std::string msg;

  // Create a log binner with 10 bins spanning the interval [1, exp(15.)):
  LogBinner binner(1., exp(15.), 10);

  // Make sure there are really 10 bins:
  msg = "LogBinner::getNumBins()";
  if (10 != binner.getNumBins()) {
    std::cerr << msg << " returned " << binner.getNumBins() << ", not 10" << std::endl;
    m_failed = true;
  }

  // Make sure values are correctly classified:
  msg = "LogBinner::computeIndex(";
  for (long ii = 0; ii < 10; ++ii) {
    double value = .9999999 * exp((ii + 1) * 15. / 10.);
    long index = binner.computeIndex(value);
    if (index != ii) {
      std::cerr << msg << value << ") returned " << index << ", not " << ii << std::endl;
      m_failed = true;
    }
  }

  // Left endpoint should be included:
  long index = binner.computeIndex(1.);
  if (0 != index) {
    std::cerr << msg << "1.) returned " << index << ", not 0" << std::endl;
    m_failed = true;
  }

  // Right endpoint should be excluded:
  index = binner.computeIndex(exp(15.));
  if (-1 != index) {
    std::cerr << msg << "exp(15.) returned " << index << ", not -1" << std::endl;
    m_failed = true;
  }

  // Left of left endpoint should return index < 0:
  index = binner.computeIndex(0.);
  if (0 <= index) {
    std::cerr << msg << "0.) returned " << index << ", which is >= 0" << std::endl;
    m_failed = true;
  }

  // Right of right endpoint should return index < 0:
  index = binner.computeIndex(1.000001 * exp(15.));
  if (0 <= index) {
    std::cerr << msg << "1.000001 * exp(15.)) returned " << index << ", which is >= 0" << std::endl;
    m_failed = true;
  }
}

void EvtBinTest::testOrderedBinner() {
  std::string msg = "OrderedBinner::OrderedBinner(...)";

  OrderedBinner::IntervalCont_t intervals;

  // Create intervals with bad ordering within a given interval.
  intervals.push_back(Binner::Interval(0., .1));
  intervals.push_back(Binner::Interval(.15, .25));
  intervals.push_back(Binner::Interval(.25, .24));
  intervals.push_back(Binner::Interval(.30, .45));

  try {
    OrderedBinner binner(intervals);
    std::cerr << msg << " did not throw when given an interval whose beginning value is > ending value" << std::endl;
    m_failed = true;
  } catch (const std::exception &) {
  }

  intervals.clear();

  // Create intervals with bad ordering between two subsequent intervals.
  intervals.push_back(Binner::Interval(0., .1));
  intervals.push_back(Binner::Interval(.15, .25));
  intervals.push_back(Binner::Interval(.24, .29));
  intervals.push_back(Binner::Interval(.30, .45));

  try {
    OrderedBinner binner(intervals);
    std::cerr << msg << " did not throw when given two sequential intervals which are not in order" << std::endl;
    m_failed = true;
  } catch (const std::exception &) {
  }

  // Finally, create a legitimate set of intervals.
  intervals.clear();

  intervals.push_back(Binner::Interval(0., .1));
  intervals.push_back(Binner::Interval(.15, .25));
  intervals.push_back(Binner::Interval(.30, .45));
  intervals.push_back(Binner::Interval(.50, .60));
  intervals.push_back(Binner::Interval(.60, .72));

  try {
    OrderedBinner binner(intervals);

    msg = "OrderedBinner::computeIndex(";

    double value;
    long index;

    // Make sure all values are classified correctly.
    // A value less than the first bin.
    value = -.01;
    index = binner.computeIndex(value);
    if (0 <= index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not a negative index" << std::endl;
    }
    
    // A value greater than the last bin.
    value = 1.;
    index = binner.computeIndex(value);
    if (0 <= index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not a negative index" << std::endl;
    }
    
    // A value on the leading edge of the first bin.
    value = 0.;
    index = binner.computeIndex(value);
    if (0 != index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not 0" << std::endl;
    }

    // A value on the trailing edge of the last bin.
    value = .72;
    index = binner.computeIndex(value);
    if (0 <= index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not a negative index" << std::endl;
    }

    // A value in the 0th bin.
    value = .05;
    index = binner.computeIndex(value);
    if (0 != index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not 0" << std::endl;
    }
    
    // A value in the 1st bin.
    value = .17;
    index = binner.computeIndex(value);
    if (1 != index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not 1" << std::endl;
    }
    
    // A value in the 2st bin.
    value = .30;
    index = binner.computeIndex(value);
    if (2 != index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not 2" << std::endl;
    }
    
    // A value in the 3rd bin.
    value = .55;
    index = binner.computeIndex(value);
    if (3 != index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not 3" << std::endl;
    }
    
    // A value in the 4th bin.
    value = .60;
    index = binner.computeIndex(value);
    if (4 != index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not 4" << std::endl;
    }
    
    // A value between bins.
    value = .25;
    index = binner.computeIndex(value);
    if (0 <= index) {
      m_failed = true;
      std::cerr << msg << value << ") returned " << index << ", not a negative index" << std::endl;
    }
  
  } catch (const std::exception &) {
    std::cerr << msg << " threw when given a set of intervals which are legal (i.e. in order)" << std::endl;
    m_failed = true;
  }
}

void EvtBinTest::testHealpixBinner() {
  std::string msg = "HealpixBinner::HealpixBinner(...)";

  //test computation of number of bins, for a given order
  for(int order=0;order!=13;order++){
    HealpixBinner binner("NESTED", order, true, "Binner");
    long nbins= binner.getNumBins();
    long int nside = pow((long double)2,order);
    long true_nbins=12*nside*nside;
    if(nbins!=true_nbins) {
      m_failed=true;
      std::cerr << msg << order << ") returned " << nbins << ", instead of " <<true_nbins<< std::endl;
    }
  }
  try { HealpixBinner binner("NESTED", -1, true, "Binner");}
  catch(const std::exception & x) {
    m_os.info() << "Expected: failed to create a Healpix Binner with order -1 : " << x.what() << std::endl;
  }
  try { HealpixBinner binner("NESTED", 13, true, "Binner");}
  catch(const std::exception & x) {
    m_os.info() << "Expected: failed to create a Healpix Binner with order 13 : " << x.what() << std::endl;
  }
}

void EvtBinTest::testHealpixMap() {
  Gti gti(m_ft1_file);
  //no energy binning case:
  LogBinner energy_binner(m_e_min, m_e_max, 0., "ENERGY");
  HealpixMap healpix_map(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA",
			 "RING", 3, false, energy_binner, energy_binner, true,gti);
  healpix_map.binInput();
  healpix_map.writeOutput("test_evtbin", "test.healmap");

  //energy binning case:
  LogBinner energy_binner2(m_e_min, m_e_max, 10, "ENERGY");
  HealpixMap healpix_cube(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA",
			 "RING", 3, true, energy_binner2, energy_binner2, true,gti);
  healpix_cube.binInput();
  healpix_cube.writeOutput("test_evtbin", "test.healcube");


}

void EvtBinTest::testHist1D() {
  std::string msg = "Hist1D";

  // Create a linear binner with bin width == 15. spanning the interval [0, 100):
  LinearBinner binner(0., 100., 15.);

  // Create a histogram using this binner:
  Hist1D lin_hist(binner);

  // Populate this histogram, starting from right of the right endpoint, going to left of left endpoint:
  for (int ii = 100; ii >= -1; --ii) lin_hist.fillBin(ii);

  // Last bin has 5 fewer values because the interval is not an integer multiple of the bin size.
  // Below it will be checked that each bin has the same number of counts, so pad it out here.
  for (int ii = 0; ii < 5; ++ii) lin_hist.fillBin(97.);

  // Check whether each bin has the right number:
  int bin_num = 0;
  for (Hist1D::ConstIterator itor = lin_hist.begin(); itor != lin_hist.end(); ++itor, ++bin_num) {
    if (15 != *itor) {
      std::cerr << msg << "'s bin number " << bin_num << " has " << *itor << " counts, not 15" << std::endl;
      m_failed = true;
    }
  }
}

void EvtBinTest::testHist2D() {
  std::string msg = "Hist2D";

  // Create a linear binner with bin width == 10 spanning the interval [0, 100):
  LinearBinner binner1(0., 100., 10.);

  // Create a log binner with 10 bins spanning the interval [1, exp(10.)):
  LogBinner binner2(1., exp(10.), 10);

  // Create a histogram using these binners:
  Hist2D hist(binner1, binner2);

  // Populate this histogram, starting from right of the right endpoint, going to left of left endpoint:
  for (int ii = 100; ii >= -1; --ii) {
    for (long jj = 0; jj < 10; ++jj) {
      double value = .9999999 * exp(double(jj + 1));
      hist.fillBin(ii, value);
    }
  }

  // Check whether each bin has the right number:
  int bin_num1 = 0;
  for (Hist2D::ConstIterator1 itor1 = hist.begin(); itor1 != hist.end(); ++itor1, ++bin_num1) {
    int bin_num2 = 0;
    for (Hist2D::ConstIterator2 itor2 = itor1->begin(); itor2 != itor1->end(); ++itor2, ++bin_num2) {
      if (10 != *itor2) {
        std::cerr << msg << "'s bin number (" << bin_num1 << ", " << bin_num2 << ") has " <<
          *itor2 << " counts, not 10" << std::endl;
        m_failed = true;
      }
    }
  }
}

void EvtBinTest::testLightCurve() {
  // Good time interval from event file.
  Gti gti(m_ft1_file);

  // Create light curve object.
  LightCurve lc(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA",
    LinearBinner(m_t_start, m_t_stop, (m_t_stop - m_t_start) * .01, "TIME"), gti);

  // Fill the light curve.
  lc.binInput();

  // Write the light curve to an output file.
  lc.writeOutput("test_evtbin", "LC1.lc");

  // Good time interval for GBM data.
  Gti gbm_gti;
  gbm_gti.insertInterval(m_gbm_t_start, m_gbm_t_stop);

  // Create light curve object for GBM data.
  LightCurve gbm_lc(m_gbm_file, "EVENTS", "", "SC_DATA",
    LinearBinner(m_gbm_t_start, m_gbm_t_stop, (m_gbm_t_stop - m_gbm_t_start) * .01, "TIME"), gbm_gti);

  // Fill the light curve.
  gbm_lc.binInput();

  // Write the light curve to an output file.
  gbm_lc.writeOutput("test_evtbin", "GBMLC1.lc");

}

void EvtBinTest::testSingleSpectrum() {
  // Test creating spectrum from LAT data.
  // Create binner used both for energy bins and for ebounds definition.
  LogBinner energy_binner(m_e_min, m_e_max, 100, "ENERGY");

  // Good time interval from event file.
  Gti gti(m_ft1_file);

  // Create spectrum object.
  SingleSpec spectrum(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA", energy_binner, energy_binner, gti);

  // Fill the spectrum.
  spectrum.binInput();

  // Write the spectrum to an output file.
  spectrum.writeOutput("test_evtbin", "PHA1.pha");

  // Test creating spectrum from GBM data.
  // Start with GBM configuration.
  std::auto_ptr<BinConfig> config(BinConfig::create(m_gbm_file));

  // Set parameters.
  st_app::AppParGroup & pars(getParGroup());
  pars["evfile"] = m_gbm_file;
  std::auto_ptr<Binner> gbm_energy_binner(config->createEnergyBinner(pars));
  std::auto_ptr<Binner> gbm_ebounds(config->createEbounds(pars));
  std::auto_ptr<Gti> gbm_gti(config->createGti(pars));

  SingleSpec gbm_spectrum(m_gbm_file, "EVENTS", "", "SC_DATA", *gbm_energy_binner, *gbm_ebounds, *gbm_gti);

  // Fill the spectrum.
  gbm_spectrum.binInput();

  // Write the spectrum to an output file.
  std::string output_file = "GBMPHA1.pha";
  gbm_spectrum.writeOutput("test_evtbin", output_file);

  // Make sure this spectrum does not contain a 2D histogram.
  try {
    gbm_spectrum.getHist2D();
    m_failed = true;
    std::cerr <<
      "Unexpected: in testSingleSpectrum, getHist2D did not throw exception when called for what should be a 1D histogram." <<
      std::endl;
  } catch (const std::exception &) {
    // OK, this is supposed to throw.
  }

  // Get the output data from the spectrum and compare it to the output file to make sure it's valid.
  const Hist1D & hist(gbm_spectrum.getHist1D());

  std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(output_file, "SPECTRUM"));
  if (0 == table->getNumRecords()) {
    m_failed = true;
    std::cerr << "Unexpected: in testSingleSpectrum, output file " << output_file << " has 0 records." << std::endl;
  } else if (hist.end() - hist.begin() != table->getNumRecords()) {
    m_failed = true;
    std::cerr << "Unexpected: in testSingleSpectrum, output file " << output_file << " has " << table->getNumRecords() <<
      " records, not " << int(hist.end() - hist.begin()) << ", as expected." << std::endl;
  } else {
    tip::Table::ConstIterator table_itor = table->begin();
    bool discrepancy = false;
    for (Hist1D::ConstIterator itor = hist.begin(); !discrepancy && itor != hist.end(); ++itor, ++table_itor) {
      if (*itor != (*table_itor)["COUNTS"].get()) discrepancy = true;
    }
    if (discrepancy) {
      m_failed = true;
      std::cerr << "Unexpected: in testSingleSpectrum, output file " << output_file << " disagrees with histogram." << std::endl;
    }
  }
}

void EvtBinTest::testMultiSpectra() {
  // Create binner used both for energy bins and for ebounds definition.
  LogBinner energy_binner(m_e_min, m_e_max, 100, "ENERGY"), LogBinner(m_e_min, m_e_max, 100, "ENERGY");

  // Good time interval from event file.
  Gti gti(m_ft1_file);

  // Create spectrum object.
  MultiSpec spectrum(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA",
    LinearBinner(m_t_start, m_t_stop, (m_t_stop - m_t_start) * .1, "TIME"), energy_binner, energy_binner, gti);

  // Fill the spectrum.
  spectrum.binInput();

  // Write the spectrum to an output file.
  spectrum.writeOutput("test_evtbin", "PHA2.pha");

  // Make sure this spectrum does not contain a 1D histogram.
  try {
    spectrum.getHist1D();
    m_failed = true;
    std::cerr <<
      "Unexpected: in testMultiSpectra, getHist1D did not throw exception when called for what should be a 1D histogram." <<
      std::endl;
  } catch (const std::exception &) {
    // OK, this is supposed to throw.
  }

  // Make sure we can extract a 2D histogram from spectrum.
  spectrum.getHist2D();
}

void EvtBinTest::testCountMap() {
  // Create map object with invalid projection.
  try {
// Disabling this test because currently astro doesn't detect the error.
//    CountMap count_map(240., 40., "GARBAGE", 400, 200, .1, 0., true, "RA", "DEC");
//    std::cerr << "CountMap's constructor did not throw an exception when given an invalid projection type" << std::endl;
//    m_failed = true;
  } catch (const std::exception &) {
    // OK, supposed to fail.
  }

  // Good time interval from event file.
  Gti gti(m_ft1_file);

  // Create count map object.
  CountMap count_map(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA", 83.4, 22.0, "AIT", 100, 100, .1, 0., false,
    "RA", "DEC", gti);

  // Fill the count map.
  count_map.binInput();

  // Write the count map to an output file.
  count_map.writeOutput("test_evtbin", "CM2.fits");
}

void EvtBinTest::testCountCube() {

  // Test creating spectrum from LAT data.
  // Create binner used both for energy bins and for ebounds definition.
  LogBinner energy_binner(m_e_min, m_e_max, 100, "ENERGY");

  // Good time interval from event file.
  Gti gti(m_ft1_file);

  // Create the count cube object.
  CountCube count_cube(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA",
		       83.4, 22.0, "AIT", 100, 100, .1, 0., false,
		       "RA", "DEC", energy_binner, energy_binner, gti);

  // Fill the count cube.
  count_cube.binInput();

  // Write the count cube to an output file.
  count_cube.writeOutput("test_evtbin", "test.ccube");

}

void EvtBinTest::testBinConfig() {
  m_os.setMethod("testBinConfig");

  // Get parameters.
  st_app::AppParGroup & par_group = getParGroup("test_evtbin");

  Binner * binner = 0;

  try {
    // Set the input file/table.
    par_group["evfile"] = m_ft1_file;
    par_group["evtable"] = "EVENTS";

    // Set name of time field to be something strange.
    par_group["tfield"] = "WackyTime";

    // First test the simple case.
    par_group["tbinalg"] = "LIN";

    // Use unlikely values for the other binning parameters.
    par_group["tstart"] = -177.;
    par_group["tstop"] = -100.;
    par_group["dtime"] = 7.;

    // Save these parameters.
    par_group.Save();

    // Try to find a binner configuration. This should fail because no prototypes were loaded yet.
    try {
      std::auto_ptr<BinConfig> config(BinConfig::create(par_group["evfile"]));
      m_failed = true;
      m_os.err() << "Unexpected: created a BinConfig before prototypes were loaded." << std::endl;
    } catch (const std::exception & x) {
      m_os.info() << "Expected: failed to create a BinConfig before prototypes were loaded: " << x.what() << std::endl;
    }

    // Load standard mission/instrument configurations.
    BinConfig::load();

    // Create a configuration object.
    std::auto_ptr<BinConfig> config(BinConfig::create(par_group["evfile"]));

    // Test prompting for time binner parameters. Since they're all hidden, their values should just be
    // the same as the values just assigned above.
    config->timeParPrompt(par_group);

    // Make sure the value set above DID take.
    if (0 != par_group["tfield"].Value().compare("WackyTime")) {
      m_failed = true;
      std::cerr << "BinConfig::timeParPrompt got name " << par_group["tfield"].Value() << ", not WackyTime" << std::endl;
    }

    // Test creating the time binner.
    binner = config->createTimeBinner(par_group);

    // Name of binner should be the value from the tfield parameter.
    if (0 != binner->getName().compare("WackyTime")) {
      m_failed = true;
      std::cerr << "BinConfig::createTimeBinner created a binner named " << binner->getName() << ", not WackyTime" << std::endl;
    }

    // The number of bins should also be consistent with the parameter file.
    if (11 != binner->getNumBins()) {
      m_failed = true;
      std::cerr << "BinConfig::createTimeBinner created a binner with " << binner->getNumBins() << " bins, not 11" << std::endl;
    }

    // Get absolute O/S dependent limit on precision.
    const double epsilon = std::numeric_limits<double>::epsilon();

    // The first bin should begin with tstart.
    if (epsilon < fabs((binner->getInterval(0).begin() + 177.) / 177.)) {
      m_failed = true;
      std::cerr << "BinConfig::createTimeBinner created a binner whose first bin begins with " <<
        binner->getInterval(0).begin() << " not -177." << std::endl;
    }

    // The last bin should end with tstop.
    if (epsilon < fabs((binner->getInterval(binner->getNumBins() - 1).end() + 100.) / 100.)) {
      m_failed = true;
      std::cerr << "BinConfig::createTimeBinner created a binner whose last bin ends with " <<
        binner->getInterval(binner->getNumBins() - 1).end() << " not -100." << std::endl;
    }

    delete binner; binner = 0;

    // Now test the logarithmic case with energy bins.
    // Set name of energy field to be something strange.
    par_group["efield"] = "WackyEnergy";

    // First test the simple case.
    par_group["ebinalg"] = "LOG";

    // Use unlikely values for the other binning parameters.
    par_group["emin"] = 1.e-7;
    par_group["emax"] = 1.;
    par_group["enumbins"] = 7;

    // Save these parameters.
    par_group.Save();

    // Prompt for energy values. Again, they're all hidden.
    config->energyParPrompt(par_group);

    // Test creating the energy binner.
    binner = config->createEnergyBinner(par_group);

    // Name of binner should be the value from the efield parameter.
    if (0 != binner->getName().compare("WackyEnergy")) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner created a binner named " << binner->getName() << ", not WackyEnergy" << std::endl;
    }

    // The number of bins should also be consistent with the parameter file.
    if (7 != binner->getNumBins()) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner created a binner with " << binner->getNumBins() << " bins, not 7" << std::endl;
    }

    // The first bin should begin with emin.
    if (epsilon < fabs((binner->getInterval(0).begin() - 1.e-7) / 1.e-7)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner created a binner whose first bin begins with " <<
        binner->getInterval(0).begin() << " not 1.e-7" << std::endl;
    }

    // The first bin should end with one order of magnitude more than emin.
    if (epsilon * 100 < fabs((binner->getInterval(0).end() - 1.e-6) / 1.e-6)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner created a binner whose first bin ends with " <<
        binner->getInterval(0).end() << " not 1.e-6" << std::endl;
    }

    // The last bin should begin with one order of magnitude less than emax.
    if (epsilon * 100 < fabs((binner->getInterval(binner->getNumBins() - 1).begin() - .1) / .1)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner created a binner whose last bin begins with " <<
        binner->getInterval(binner->getNumBins() - 1).begin() << " not .1" << std::endl;
    }

    // The last bin should end with emax.
    if (epsilon < fabs((binner->getInterval(binner->getNumBins() - 1).end() - 1.) / 1.)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner created a binner whose last bin ends with " <<
        binner->getInterval(binner->getNumBins() - 1).end() << " not 1." << std::endl;
    }

    delete binner; binner = 0;

    // Now test the bin file case with energy bins.
    par_group["ebinalg"] = "FILE";
    par_group["ebinfile"] = facilities::commonUtilities::joinPath(m_data_dir, "energybins.fits");

    // Save these parameters.
    par_group.Save();

    // Prompt for energy values. Again, they're all hidden.
    config->energyParPrompt(par_group);

    // Test creating the energy binner.
    binner = config->createEnergyBinner(par_group);

    // The number of bins should also be consistent with the bin definition file.
    if (1024 != binner->getNumBins()) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner created a binner with " << binner->getNumBins() << " bins, not 1024" << std::endl;
    }

    // The beginning of the first bin should match the file contents.
    if (epsilon < fabs((binner->getInterval(0).begin() - 30.) / 30.)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner with bin def file created a binner whose first bin begins with " <<
        binner->getInterval(0).begin() << " not 30." << std::endl;
    }

    float float_epsilon = std::numeric_limits<float>::epsilon();
    // The end value of the first bin should match the file contents.
    if (float_epsilon < fabs((binner->getInterval(0).end() - 30.20306) / 30.20306)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner with bin def file created a binner whose first bin ends with " <<
        binner->getInterval(0).end() << " not 30.20306" << std::endl;
    }

    // The beginning value of the last bin should match the file contents.
    if (float_epsilon < fabs((binner->getInterval(binner->getNumBins() - 1).begin() - 29798.306) / 29798.306)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner with bin def file created a binner whose last bin begins with " <<
        binner->getInterval(binner->getNumBins() - 1).begin() << " not 29798.306" << std::endl;
    }

    // The end value of the last bin should match the file contents.
    if (float_epsilon < fabs((binner->getInterval(binner->getNumBins() - 1).end() - 30000.) / 30000.)) {
      m_failed = true;
      std::cerr << "BinConfig::createEnergyBinner with bin def file created a binner whose last bin ends with " <<
        binner->getInterval(binner->getNumBins() - 1).end() << " not 30000." << std::endl;
    }

  } catch (const std::exception & x) {
    m_failed = true;
    std::cerr << "testBinConfig encountered an unexpected error: " << x.what() << std::endl;
  }

  delete binner;

}

void EvtBinTest::testGti() {
  Gti gti1;
  gti1.insertInterval(1., 2.);

  Gti gti2;
  gti2.insertInterval(2., 3.);

  Gti result = gti1 & gti2;
  if (result.begin() != result.end()) {
    std::cerr << "testGti found gti1 overlaps gti2 but they should be disjoint" << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti2 & gti1;
  if (result.begin() != result.end()) {
    std::cerr << "testGti found gti2 overlaps gti1 but they should be disjoint" << std::endl;
    m_failed = true;
  }

  Gti gti3;
  gti3.insertInterval(1.5, 1.75);

  result = Gti();
  result = gti1 & gti3;
  if (result.begin() == result.end()) {
    std::cerr << "testGti found gti1 does not overlap gti3 but they should overlap" << std::endl;
    m_failed = true;
  } else if (result.begin()->first != 1.5 || result.begin()->second != 1.75) {
    std::cerr << "testGti found gti1 & gti3 == [" << result.begin()->first << ", " << result.begin()->second <<
      "], not [1.5, 1.75]" << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti3 & gti1;
  if (result.begin() == result.end()) {
    std::cerr << "testGti found gti3 does not overlap gti1 but they should overlap" << std::endl;
    m_failed = true;
  } else if (result.begin()->first != 1.5 || result.begin()->second != 1.75) {
    std::cerr << "testGti found gti3 & gti1 == [" << result.begin()->first << ", " << result.begin()->second <<
      "], not [1.5, 1.75]" << std::endl;
    m_failed = true;
  }

  Gti gti4;
  gti4.insertInterval(1.5, 2.5);

  result = Gti();
  result = gti1 & gti4;
  if (result.begin() == result.end()) {
    std::cerr << "testGti found gti1 does not overlap gti4 but they should overlap" << std::endl;
    m_failed = true;
  } else if (result.begin()->first != 1.5 || result.begin()->second != 2.0) {
    std::cerr << "testGti found gti1 & gti4 == [" << result.begin()->first << ", " << result.begin()->second <<
    "], not [1.5, 2.0]" << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti4 & gti1;
  if (result.begin() == result.end()) {
    std::cerr << "testGti found gti4 does not overlap gti1 but they should overlap" << std::endl;
    m_failed = true;
  } else if (result.begin()->first != 1.5 || result.begin()->second != 2.0) {
    std::cerr << "testGti found gti4 & gti1 == [" << result.begin()->first << ", " << result.begin()->second <<
      "], not [1.5, 2.0]" << std::endl;
    m_failed = true;
  }

  // Now two GTIs with multiple entries.
  Gti gti5;
  gti5.insertInterval(1., 2.);
  gti5.insertInterval(3., 4.);
  gti5.insertInterval(5., 6.);
  gti5.insertInterval(9., 10.);

  Gti gti6;
  gti6.insertInterval(2.5, 3.5);
  gti6.insertInterval(3.75, 5.1);
  gti6.insertInterval(5.3, 5.5);
  gti6.insertInterval(6.5, 7.5);
  gti6.insertInterval(8.5, 9.5);
  gti6.insertInterval(10.5, 11.5);

  Gti correct_result;
  correct_result.insertInterval(3., 3.5);
  correct_result.insertInterval(3.75, 4.);
  correct_result.insertInterval(5., 5.1);
  correct_result.insertInterval(5.3, 5.5);
  correct_result.insertInterval(9., 9.5);

  result = Gti();
  result = gti5 & gti6;
  if (result != correct_result) {
    std::cerr << "testGti found gti5 & gti6 did not return expected result" << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti6 & gti5;
  if (result != correct_result) {
    std::cerr << "testGti found gti6 & gti5 did not return expected result" << std::endl;
    m_failed = true;
  }

  // Check ONTIME computation.
  double on_time = correct_result.computeOntime();
  double expected_on_time = 1.55;
  const double tolerance = 1.e-12;
  if (tolerance < fabs(expected_on_time - on_time)) {
    std::cerr << "testGti: computeOntime returned " << on_time << ", not " << expected_on_time << " as expected" << std::endl;
    m_failed = true;
  }

  Gti gti(m_ft1_file);

  // Create light curve object.
  LightCurve lc(m_ft1_file, "EVENTS", m_ft2_file, "SC_DATA",
    LinearBinner(m_t_start, m_t_stop, (m_t_stop - m_t_start) * .01, "TIME"), gti);

  // Get absolute O/S dependent limit on precision.
  //  const double epsilon = std::numeric_limits<double>::epsilon();
  const double epsilon = 1.0e-10;

  const Gti & lc_gti = lc.getGti();
  const int n_gti = 13;
  if (n_gti != lc_gti.getNumIntervals()) {
    std::cerr << "testGti: ERROR: read GTI from test ft1 file with " << lc_gti.getNumIntervals() << ", not " << n_gti << " GTI" << std::endl;
    m_failed = true;
  } else {
    // The first interval should agree exactly with the start interval read from the file.
    const double gti1_duration = 4399.12953472137;
    if (m_t_start != lc_gti.begin()->first || epsilon < fabs((m_t_start + gti1_duration - lc_gti.begin()->second) / (m_t_start + gti1_duration))) {
      std::cerr << "testGti: ERROR: read GTI from test ft1 file with values [" << lc_gti.begin()->first << ", " <<
        lc_gti.begin()->second << ", not [" << m_t_start << ", " << m_t_start + gti1_duration << "], as expected." << std::endl;
      m_failed = true;
    }

    Gti::ConstIterator last = lc_gti.end();
    --last;
    // The last interval should agree less exactly with the stop interval read from the file, because the gti was recut.
    const double gti2_duration = 1.692876845e+3;
    if (m_t_stop != last->second || 10. * epsilon < fabs((m_t_stop - gti2_duration - last->first) / (m_t_stop - gti2_duration))) {
      std::cerr << "testGti: ERROR: read GTI from test ft1 file with values [" << last->first << ", " << last->second <<
        ", not [" << m_t_stop - gti2_duration << ", " << m_t_stop << "], as expected." << std::endl;
      m_failed = true;
    }
  }

  // Check ONTIME computation from light curve.
  expected_on_time = 84307.6959530115;
  on_time = lc_gti.computeOntime();
  // Error in ontime is compounded by the large number of bins used in the light curve.
  if (1000. * epsilon < fabs((expected_on_time - on_time) / expected_on_time)) {
    std::cerr << "testGti: computeOntime returned " << on_time << ", not " << expected_on_time << " as expected" << std::endl;
    m_failed = true;
  }

  Gti spec_gti("PHA1.pha", "GTI");

  Gti::ConstIterator gti_pos = spec_gti.begin();

  const double gti_start = 223850000.0;
  const double gti_stop = 223854399.129535;

  std::cerr.precision(24);

  // Interval == GTI.
  double fract = spec_gti.getFraction(gti_start, gti_stop, gti_pos);
  if (fabs(fract - 1.0) > epsilon) {
    std::cerr << "Unexpected: testGti: Interval == GTI, getFraction returned " << fract << ", not 1" << std::endl;
    m_failed = true;
//  if (gti_pos !=spec_gti.end())
//    std::cerr << "Unexpected: testGti: Interval == GTI, iterator was not incremented" << std::endl;
  }

  // Interval < GTI.
  gti_pos = spec_gti.begin();
  fract = spec_gti.getFraction(gti_start - 1.e3, gti_start, gti_pos);
  if (fract != 0.) {
    std::cerr << "Unexpected: testGti: Interval < GTI, getFraction returned " << fract << ", not 0" << std::endl;
    m_failed = true;
  }
  if (gti_pos != spec_gti.begin()) {
    std::cerr << "Unexpected: testGti: Interval < GTI, iterator was incremented" << std::endl;
    m_failed = true;
  }

  // Interval > GTI.
  gti_pos = spec_gti.begin();
  Gti::ConstIterator gti_pos2 = gti_pos;
  gti_pos2++;
  fract = spec_gti.getFraction(gti_stop + .0001, gti_stop + 1.e3, gti_pos);
  if (fract != 0.) {
    std::cerr << "Unexpected: testGti: Interval > GTI, getFraction returned " << fract << ", not 0" << std::endl;
    m_failed = true;
  }



  //for (int i = 0; i < (n_gti - 1); i++) { ++gti_pos; }
  //  if (++gti_pos != spec_gti.end()) {
  if (gti_pos != gti_pos2) {
    std::cerr << "Unexpected: testGti: Interval > GTI, iterator was not incremented" << std::endl;
    m_failed = true;
  }

  // Interval starts before GTI, goes halfway through.
  gti_pos = spec_gti.begin();
  double gti_width = gti_stop - gti_start;
  fract = spec_gti.getFraction(gti_start - gti_width * .5, gti_stop - gti_width * .5, gti_pos);
  if (epsilon < fabs((fract - .5) / .5)) {
    std::cerr << "Unexpected: testGti: Interval starts before GTI, getFraction returned " << fract << ", not .5" << std::endl;
    m_failed = true;
  }
  if (gti_pos != spec_gti.begin()) {
    std::cerr << "Unexpected: testGti: Interval starts before GTI, iterator was incremented" << std::endl;
    m_failed = true;
  }

  // Interval starts ~10%-way through GTI, goes past end.
  gti_pos = spec_gti.begin();
  const double fw = 0.1;
  const double omfw = 1.0 - fw;
  fract = spec_gti.getFraction(gti_start + gti_width * fw,
                               gti_stop + gti_width * fw, gti_pos);
  if (epsilon < fabs((fract - omfw) / omfw)) {
    std::cerr << "Unexpected: testGti: GTI starts before interval, "
      "getFraction returned " << fract << ", not " << omfw << std::endl;
    m_failed = true;
  }




  if (gti_pos != gti_pos2) {
    std::cerr << "Unexpected: testGti: GTI starts before interval, iterator was not incremented" << std::endl;
    m_failed = true;
  }

  // Interval contained within GTI.
  gti_pos = spec_gti.begin();
  fract = spec_gti.getFraction(gti_start + 100., gti_stop - 100., gti_pos);
  if (fract != 1.) {
    std::cerr << "Unexpected: testGti: Interval contained within GTI, getFraction returned " << fract << ", not 1." << std::endl;
    m_failed = true;
  }
  if (gti_pos != spec_gti.begin()) {
    std::cerr << "Unexpected: testGti: Interval contained within GTI, iterator was incremented" << std::endl;
    m_failed = true;
  }

  // GTI contained within interval.
  gti_pos = spec_gti.begin();
  fract = spec_gti.getFraction(gti_start - gti_width * .5, gti_stop + gti_width * .1, gti_pos);
  if (epsilon < fabs((fract -.625) / .625)) {
    std::cerr << "Unexpected: testGti: GTI contained within interval, getFraction returned " << fract << ", not .625" << std::endl;
    m_failed = true;
  }
  if (gti_pos != gti_pos2) {
    std::cerr << "Unexpected: testGti: GTI contained within interval, iterator was not incremented" << std::endl;
    m_failed = true;
  }

  // Test logical or of two gtis.
  gti5 = Gti();
  gti5.insertInterval(1., 2.);
  gti5.insertInterval(3., 4.);
  gti5.insertInterval(5., 6.);

  gti6 = Gti();
  gti6.insertInterval(2., 2.5);
  gti6.insertInterval(3.5, 4.5);
  gti6.insertInterval(6.5, 7.5);

  correct_result = Gti();
  correct_result.insertInterval(1., 2.5);
  correct_result.insertInterval(3., 4.5);
  correct_result.insertInterval(5., 6.);
  correct_result.insertInterval(6.5, 7.5);

  // Try or-ing two interleaving Gtis.
  result = Gti();
  result = gti5 | gti6;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after gti5 | gti6, result was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti6 | gti5;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after gti6 | gti5, result was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  // Try or-ing two empty Gtis
  gti3 = Gti();
  gti4 = Gti();
  correct_result = Gti();
  result = Gti();
  result = gti3 | gti4;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after gti3 | gti4, result was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  // Try or-ing one empty, one non-empty Gti.
  correct_result = gti6;
  result = Gti();
  result = gti4 | gti6;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after gti4 | gti6, result was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  // Try or-ing one non-empty, one empty Gti.
  correct_result = gti5;
  result = Gti();
  result = gti5 | gti3;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after gti5 | gti3, result was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  gti3 = Gti();
  gti3.insertInterval(10., 15.);
  gti4 = Gti();
  gti4.insertInterval(10., 20.);
  correct_result = gti4;
  result = Gti();
  result = gti3 | gti4;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting both gti to start at the same time, result of gti3 | gti4 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti4 | gti3;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting both gti to start at the same time, result of gti4 | gti3 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  gti3 = Gti();
  gti3.insertInterval(15., 20.);
  gti4 = Gti();
  gti4.insertInterval(10., 20.);

  correct_result = gti4;
  result = Gti();
  result = gti3 | gti4;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting both gti to stop at the same time, result of gti3 | gti4 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti4 | gti3;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting both gti to stop at the same time, result of gti4 | gti3 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  gti3 = Gti();
  gti3.insertInterval(15., 18.);
  gti4 = Gti();
  gti4.insertInterval(10., 20.);

  correct_result = gti4;
  result = Gti();
  result = gti3 | gti4;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting gti3 to be contained completely inside gti4, result of gti3 | gti4 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti4 | gti3;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting gti3 to be contained completely inside gti4, result of gti4 | gti3 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  correct_result = gti3;
  result = Gti();
  result = gti3 & gti4;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting gti3 to be contained completely inside gti4, result of gti3 & gti4 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

  result = Gti();
  result = gti4 & gti3;
  if (result != correct_result) {
    std::cerr << "Unexpected: testGti: after setting gti3 to be contained completely inside gti4, result of gti4 & gti3 was:\n" <<
      result << "\n, not\n" << correct_result << "\nas expected." << std::endl;
    m_failed = true;
  }

}

void EvtBinTest::testConstSnBinner() {
  m_os.setMethod("testConstSnBinner()");

  // Create test binner object.
  ConstSnBinner binner(1., 101., 5., 1., 25.);

  // Check value before first bin.
  long index = binner.computeIndex(0.);
  if (-1 != index) {
    m_failed = true;
    m_os.err() << "In first binner test, index of 0. was " << index << ", not -1" << std::endl;
  }

  // Check value exactly at the first bin.
  index = binner.computeIndex(1.);
  if (-1 != index) {
    m_failed = true;
    m_os.err() << "In first binner test, index of 1. was " << index << ", not -1" << std::endl;
  }

  // Check values up through and including the end of the last bin.
  for (double ev_time = 2.; ev_time <= 101.; ++ev_time) {
    index = binner.computeIndex(ev_time);
    long correct_index = long((ev_time - 2.)/25.);
    if (index != correct_index) {
      m_failed = true;
      m_os.err() << "In first binner test, index of " << ev_time << " was " << index << ", not " << correct_index << std::endl;
    }
  }

  // Check value after last bin.
  index = binner.computeIndex(102.);
  if (-1 != index) {
    m_failed = true;
    m_os.err() << "In first binner test, index of 102. was " << index << ", not -1" << std::endl;
  }

  // Reset the binner.
  binner = ConstSnBinner(1., 101., 5., 1., 25.);

  // Loop over some more times with larger step between each one.
  for (double ev_time = 2.; ev_time <= 101.; ev_time += 2.) {
    index = binner.computeIndex(ev_time);
    long correct_index = long((ev_time - 2.)/50.);
    if (index != correct_index) {
      m_failed = true;
      m_os.err() << "In second binner test, index of " << ev_time << " was " << index << ", not " << correct_index << std::endl;
    }
  }

  // Check value after last bin.
  index = binner.computeIndex(102.);
  if (-1 != index) {
    m_failed = true;
    m_os.err() << "In second binner test, index of 102. was " << index << ", not -1" << std::endl;
  }

  // Reset the binner.
  binner = ConstSnBinner(1., 101., 5., 1., 25.);

  // Loop over some times which are not uniformly distributed. Note upper cutoff is different from earlier tests
  // because adding the sine to the event time makes the last times later than the binner's upper bound.
  for (double ev_time = 2.; ev_time < 101.; ++ev_time) {
    double true_ev_time = ev_time + .1 * sin(ev_time);
    index = binner.computeIndex(true_ev_time);
    long correct_index = long((ev_time - 2.)/25.);
    if (index != correct_index) {
      m_failed = true;
      m_os.err() << "In third binner test, index of " << true_ev_time << " was " << index << ", not " << correct_index << std::endl;
    }
  }

  // Reset the binner, this time with a uniform background.
  std::vector<double> background_coeffs(1, 1.);
  binner = ConstSnBinner(1., 101., 5., 1., 25., background_coeffs);

  // Loop over twice as many times, but with the background of 1, the results will be the same as the second test.

  // Check value before first bin.
  index = binner.computeIndex(0.);
  if (-1 != index) {
    m_failed = true;
    m_os.err() << "In fourth binner test, index of 0. was " << index << ", not -1" << std::endl;
  }

  // Loop over some more times with smaller step between each one.
  for (double ev_time = 1.5; ev_time <= 101.; ev_time += .5) {
    index = binner.computeIndex(ev_time);
    long correct_index = long((ev_time - 1.5)/50.);
    if (index != correct_index) {
      m_failed = true;
      m_os.err() << "In fourth binner test, index of " << ev_time << " was " << index << ", not " << correct_index << std::endl;
    }
  }

  // Check value after last bin.
  index = binner.computeIndex(102.);
  if (-1 != index) {
    m_failed = true;
    m_os.err() << "In fourth binner test, index of 102. was " << index << ", not -1" << std::endl;
  }

}

void EvtBinTest::testBayesianBinner() {
  m_os.setMethod("testBayesianBinner()");

  // Create uniform bins, and highly artificial data which obviously hovers around constant values.
  OrderedBinner::IntervalCont_t intervals;
  std::vector<double> cell_pop(100);
  double change_points[] = { 0., 20., 21., 22., 27., 40., 60., 100. };
  double plateau[] = { 5., 25., 150., 260., 120., 58., 12., 0. };
  int cp = 1;
  for (int ii = 0; ii < 100; ++ii) {
    intervals.push_back(Binner::Interval(ii, ii + 1));
    if (ii >= change_points[cp]) ++cp;
    cell_pop[ii] = plateau[cp] + 3. * rand() / RAND_MAX;
  }

  // Use BayesianBinner to determine blocks, and compare with the correct answer.
  BayesianBinner binner(intervals, cell_pop.begin(), "", 6.);
  int num_bins = (sizeof(change_points) / sizeof(double)) - 1;
  if (binner.getNumBins() != num_bins) {
    m_failed = true;
    m_os.err() << "Number of Bayesian blocks found was " << binner.getNumBins() << ", not " << num_bins <<
      ", as expected." << std::endl;
    for (int ii = 0; ii != binner.getNumBins(); ++ii) {
      Binner::Interval interval(binner.getInterval(ii));
      m_os.err() << "Interval[" << ii << "] is [" << interval.begin() << ", " << interval.end() << "]" << std::endl;
    }
  } else {
    for (int ii = 0; ii != num_bins; ++ii) {
      Binner::Interval interval(binner.getInterval(ii));
      if (interval.begin() != change_points[ii] || interval.end() != change_points[ii + 1]) {
        m_failed = true;
        m_os.err() << "Interval[" << ii << "] is [" << interval.begin() << ", " << interval.end() << "], not [" <<
          change_points[ii] << ", " << change_points[ii + 1] << "], as expected." << std::endl;
      }
    }
  }
}

void EvtBinTest::testMultipleFiles() {
  using namespace st_facilities;
  m_os.setMethod("testMultipleFiles()");

  // Names of unsplit files containing events and sc data.
  std::string ev_file = m_ft1_file;
  std::string sc_file = m_ft2_file;

  // Get Gti from unsplit file.
  Gti separate_gti(ev_file);

  // Create binner used both for energy bins and for ebounds definition.
  LogBinner energy_binner(30., 200000., 100, "ENERGY");

  // Create spectrum from unsplit files.
  SingleSpec separate(ev_file, "EVENTS", sc_file, "SC_DATA", energy_binner, energy_binner, separate_gti);

  // Fill the spectrum.
  separate.binInput();

  // Write the spectrum to an output file.
  separate.writeOutput("test_evtbin", "separate_spectrum.pha");

  // Names of files containing lists of equivalent split files containing events and sc data.
  std::string ev_list_file = "@" + facilities::commonUtilities::joinPath(m_data_dir, "ft1filelist");
  std::string sc_list_file = "@" + facilities::commonUtilities::joinPath(m_data_dir, "ft2filelist");

  // Expand the contents of the list file.
  FileSys::FileNameCont tmp_input_file = FileSys::expandFileList(ev_list_file);

  // Make a vector of test files (just in case FileNameCont ever changes type.
  std::vector<std::string> input_file(tmp_input_file.begin(), tmp_input_file.end());

  // Number of events in each file.
  tip::Index_t expected_num_events = 0;
  std::vector<tip::Index_t> num_rec(input_file.size());

  // Get information about input files for verification purposes.
  std::vector<std::string>::size_type index = 0;
  for (FileSys::FileNameCont::iterator itor = input_file.begin(); itor != input_file.end(); ++itor, ++index) {
    std::auto_ptr<const tip::Table> table(tip::IFileSvc::instance().readTable(*itor, "EVENTS"));
    num_rec[index] = table->getNumRecords();
    expected_num_events += num_rec[index];
  }

  // Get the good time intervals from event file list.
  Gti merged_gti(ev_list_file);
  Gti::ConstIterator merged_pos = merged_gti.begin();

  // Confirm the gtis were merged correctly.
  Gti expected_gti;
  const int n_gti = 13;
  const double gti_start[n_gti] = {223850000.0, 223856049.158012, 223862084.197088,
                                   223868030.092777, 223873894.124928,
                                   223879825.027392, 223885778.155413,
                                   223891726.103164, 223897231.049686,
                                   223930038.17317, 223936206.009276,
                                   223942287.099427, 223948307.123155};
  const double gti_stop[n_gti] = {223854399.129535, 223860305.006878,
                                  223866272.158332, 223872392.133007,
                                  223878540.140182, 223884651.061104,
                                  223890762.015572, 223896970.086514,
                                  223929670.094502, 223934706.009866,
                                  223940596.057341,223946500.116608,
                                  223950000.0};
  for (int i = 0; i < n_gti; i++) {
    expected_gti.insertInterval(gti_start[i], gti_stop[i]);
  }
  Gti::ConstIterator expected_pos = expected_gti.begin();

  //  if (merged_gti != expected_gti) {
  const double epsilon = 1.0e-10;
  for (; merged_pos != merged_gti.end(); ++merged_pos, ++expected_pos) {
    if (fabs((merged_pos->first - expected_pos->first) / expected_pos->first) > epsilon) {
      m_failed = true;
      std::cerr << "Unexpected: testMultipleFiles: Gti first computed from merged list of event files was:\n" <<
        merged_pos->first << "\nnot:\n" << expected_pos->first << "\n, as expected." << std::endl;
    }
    if (fabs((merged_pos->second - expected_pos->second) / expected_pos->second) > epsilon) {
      m_failed = true;
      std::cerr << "Unexpected: testMultipleFiles: Gti second computed from merged list of event files was:\n" <<
        merged_pos->second << "\nnot:\n" << expected_pos->second << "\n, as expected." << std::endl;
    }
  }

  // Bin up the input to make a spectrum, and require the total number of binned counts to agree with the inputs.
  SingleSpec merged(ev_list_file, "EVENTS", sc_list_file, "SC_DATA", energy_binner, energy_binner, merged_gti);

  // Fill the spectrum.
  merged.binInput();

  // Write the spectrum to an output file.
  merged.writeOutput("test_evtbin", "merged_spectrum.pha");

  // Read the output and confirm it has the requisite properties.
  std::auto_ptr<const tip::Table> spec_table(tip::IFileSvc::instance().readTable("merged_spectrum.pha", "SPECTRUM"));
  
  double num_events = 0.;
  for (tip::Table::ConstIterator itor = spec_table->begin(); itor != spec_table->end(); ++itor) {
    num_events += (*itor)["COUNTS"].get();
  }

  if (num_events != double(expected_num_events)) {
    m_failed = true;
    m_os.err() << "Created spectrum merged_spectrum.pha has " << num_events << " events, not " << expected_num_events <<
      ", as expected." << std::endl;
  }
}

/// \brief Create factory singleton object which will create the application:
st_app::StAppFactory<EvtBinTest> g_app_factory("test_evtbin");
