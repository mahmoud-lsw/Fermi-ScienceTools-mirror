/**
 * @file test.cxx
 * @brief Tests program for Cuts class.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/dataSubselector/src/test/test.cxx,v 1.39.2.1 2015/04/23 20:50:49 jchiang Exp $
 */ 

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cstdio>

#include <fstream>
#include <stdexcept>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "st_facilities/Environment.h"
#include "st_facilities/FitsUtil.h"
#include "st_facilities/Util.h"

#include "facilities/commonUtilities.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "dataSubselector/BitMaskCut.h"
#include "dataSubselector/Cuts.h"
#include "dataSubselector/Gti.h"
#include "dataSubselector/VersionCut.h"

class DssTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(DssTests);

   CPPUNIT_TEST(test_accept2);
   CPPUNIT_TEST(compareGtis);
   CPPUNIT_TEST(updateGti);
   CPPUNIT_TEST(combineGtis);
   CPPUNIT_TEST(compareCuts);
   CPPUNIT_TEST(compareCutsWithoutGtis);
   CPPUNIT_TEST(mergeGtis);
   CPPUNIT_TEST(compareUnorderedCuts);
   CPPUNIT_TEST(cutsConstructor);
   CPPUNIT_TEST(test_SkyCone);
   CPPUNIT_TEST(test_DssFormatting);
   CPPUNIT_TEST(test_removeRangeCuts);
   CPPUNIT_TEST(test_mergeRangeCuts);
   CPPUNIT_TEST(test_BitMaskCut);
   CPPUNIT_TEST(test_VersionCut);
   CPPUNIT_TEST(test_irfName);
   CPPUNIT_TEST(test_rangeCut);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();
   void test_accept2();
   void compareGtis();
   void updateGti();
   void combineGtis();
   void compareCuts();
   void compareCutsWithoutGtis();
   void mergeGtis();
   void compareUnorderedCuts();
   void cutsConstructor();
   void test_SkyCone();
   void test_DssFormatting();
   void test_removeRangeCuts();
   void test_mergeRangeCuts();
   void test_BitMaskCut();
   void test_VersionCut();
   void test_irfName();
   void test_rangeCut();

private:

   std::string m_infile;
   std::string m_outfile;
   std::string m_outfile2;
   std::string m_evtable;

   const tip::Table * m_inputTable;
   tip::Table * m_outputTable;
   tip::Table * m_outputTable2;

   tip::Table::ConstIterator m_inputIt;
   tip::Table::Iterator m_outputIt;
   tip::Table::Iterator m_output2It;

};

#define ASSERT_EQUALS(X, Y) CPPUNIT_ASSERT(std::fabs( (X - Y)/Y ) < 1e-4)

void DssTests::setUp() {
   std::string root_path 
      = st_facilities::Environment::dataPath("dataSubselector");
   m_infile = "input_events.fits";
   m_evtable = "EVENTS";
   if (root_path != "") {
      m_infile = facilities::commonUtilities::joinPath(root_path, m_infile);
   } else {
      throw std::runtime_error("Unable to determine dataSubselector's "
                               "data path");
   }
   m_outfile = "filtered_events.fits";
   m_outfile2 = "filtered_events_2.fits";
}

void DssTests::tearDown() {
}

void DssTests::test_accept2() {
   dataSubselector::Gti gti;
   for (size_t i(0); i < 100; i++) {
      double tmin = static_cast<double>(i);
      double tmax = tmin + 0.1;
      gti.insertInterval(tmin, tmax);
   }

   CPPUNIT_ASSERT(!gti.accept2(-0.05));

   CPPUNIT_ASSERT(gti.accept2(0.05));
   CPPUNIT_ASSERT(!gti.accept2(0.15));

   CPPUNIT_ASSERT(gti.accept2(50.05));
   CPPUNIT_ASSERT(!gti.accept2(50.15));

   CPPUNIT_ASSERT(!gti.accept(-0.05));

   CPPUNIT_ASSERT(gti.accept(0.05));
   CPPUNIT_ASSERT(!gti.accept(0.15));

   CPPUNIT_ASSERT(gti.accept(50.05));
   CPPUNIT_ASSERT(!gti.accept(50.15));

   CPPUNIT_ASSERT(gti.minValue() == 0);
   CPPUNIT_ASSERT(gti.maxValue() == 99.1);
}

void DssTests::compareGtis() {
   const tip::Table * gtiTable = 
      tip::IFileSvc::instance().readTable(m_infile, "GTI");

   dataSubselector::Gti gti1(*gtiTable);
   dataSubselector::Gti gti2(m_infile);

   CPPUNIT_ASSERT(!(gti1 != gti2));

   gti1.insertInterval(100000., 100010.);

   CPPUNIT_ASSERT(gti1 != gti2);
}

void DssTests::updateGti() {
   dataSubselector::Gti gti;
   gti.insertInterval(0, 1000.);
   gti.insertInterval(1500, 2000.);
   dataSubselector::Gti new_gti = gti.applyTimeRangeCut(500., 1750.);

   double expected_values[2][2] = {{500, 1000}, {1500, 1750}};
   evtbin::Gti::ConstIterator interval;
   int i(0);
   for (interval = new_gti.begin();
        interval != new_gti.end(); ++interval, i++) {
      ASSERT_EQUALS(interval->first, expected_values[i][0]);
      ASSERT_EQUALS(interval->second, expected_values[i][1]);
   }
}

void DssTests::cutsConstructor() {
   dataSubselector::Cuts my_cuts(m_infile, m_evtable);

   CPPUNIT_ASSERT(my_cuts.size() == 2);

   std::map<std::string, double> params;
   params["ENERGY"] = 20.;
   CPPUNIT_ASSERT(!my_cuts.accept(params));
   params["ENERGY"] = 30.;
   CPPUNIT_ASSERT(!my_cuts.accept(params));
   params["ENERGY"] = 30.0001;
   CPPUNIT_ASSERT(my_cuts.accept(params));
   params["ENERGY"] = 100.;
   CPPUNIT_ASSERT(my_cuts.accept(params));
   params["ENERGY"] = 2e5;
   CPPUNIT_ASSERT(my_cuts.accept(params));
   params["ENERGY"] = 2.1e5;
   CPPUNIT_ASSERT(!my_cuts.accept(params));

   params["ENERGY"] = 100.;
   params["TIME"] = 100.;
   CPPUNIT_ASSERT(my_cuts.accept(params));
   params["TIME"] = 9e4;
   CPPUNIT_ASSERT(!my_cuts.accept(params));
}

void DssTests::test_SkyCone() {
   dataSubselector::Cuts my_cuts(m_outfile, m_evtable);
   
   std::map<std::string, double> params;

   params["RA"] = 85;
   params["DEC"] = 22;
   CPPUNIT_ASSERT(my_cuts.accept(params));

   params["DEC"] = -40;
   CPPUNIT_ASSERT(!my_cuts.accept(params));
}

void DssTests::combineGtis() {
   dataSubselector::Gti gti1, gti2;

   gti1.insertInterval(200., 300.);
   gti2.insertInterval(500., 700.);

   dataSubselector::Gti new_gti, test_gti;
   test_gti.insertInterval(200., 300.);
   test_gti.insertInterval(500., 700.);

   new_gti = gti1 | gti2;
   CPPUNIT_ASSERT(new_gti.getNumIntervals() == 2);
   CPPUNIT_ASSERT(!(new_gti != test_gti));

   new_gti = gti2 | gti1;
   CPPUNIT_ASSERT(new_gti.getNumIntervals() == 2);
   CPPUNIT_ASSERT(!(new_gti != test_gti));

   dataSubselector::Gti gti3, gti4, test_gti34;
   gti3.insertInterval(200., 500.);
   gti4.insertInterval(300., 700.);
   test_gti34.insertInterval(200., 700.);

   new_gti = gti3 | gti4;
   CPPUNIT_ASSERT(new_gti.getNumIntervals() == 1);
   CPPUNIT_ASSERT(!(new_gti != test_gti34));
   
   new_gti = gti4 | gti3;
   CPPUNIT_ASSERT(new_gti.getNumIntervals() == 1);
   CPPUNIT_ASSERT(!(new_gti != test_gti34));

   dataSubselector::Gti gti5, gti6;
   gti5.insertInterval(200., 700.);
   gti6.insertInterval(300., 500.);
   
   new_gti = gti5 | gti6;
   CPPUNIT_ASSERT(new_gti.getNumIntervals() == 1);
   CPPUNIT_ASSERT(!(new_gti != test_gti34));
   
   new_gti = gti6 | gti5;
   CPPUNIT_ASSERT(new_gti.getNumIntervals() == 1);
   CPPUNIT_ASSERT(!(new_gti != test_gti34));

   gti1.insertInterval(500, 600);
   gti1.insertInterval(900, 1000);

   gti2.insertInterval(800, 850);

   new_gti = gti1 | gti2;

   dataSubselector::Gti test_gti7;
   test_gti7.insertInterval(200, 300);
   test_gti7.insertInterval(500, 700);
   test_gti7.insertInterval(800, 850);
   test_gti7.insertInterval(900, 1000);

   CPPUNIT_ASSERT(!(new_gti != test_gti7));

   new_gti = gti2 | gti1;
   CPPUNIT_ASSERT(!(new_gti != test_gti7));

   test_gti7.insertInterval(1100, 1111);
   CPPUNIT_ASSERT(new_gti != test_gti7);
}

void DssTests::compareCuts() {
   std::string extension("EVENTS");

   if (st_facilities::Util::fileExists(m_outfile)) {
      std::remove(m_outfile.c_str());
   }
   if (st_facilities::Util::fileExists(m_outfile2)) {
      std::remove(m_outfile2.c_str());
   }
   
   tip::IFileSvc::instance().createFile(m_outfile, m_infile);
   tip::IFileSvc::instance().createFile(m_outfile2, m_infile);
   
   m_inputTable = tip::IFileSvc::instance().readTable(m_infile, extension);
   
   m_outputTable = tip::IFileSvc::instance().editTable(m_outfile, extension);
   m_outputTable->setNumRecords(m_inputTable->getNumRecords());

   m_outputTable2 = tip::IFileSvc::instance().editTable(m_outfile2, extension);
   m_outputTable2->setNumRecords(m_inputTable->getNumRecords());
   
   m_inputIt = m_inputTable->begin();
   m_outputIt = m_outputTable->begin();
   m_output2It = m_outputTable2->begin();

   tip::ConstTableRecord & input = *m_inputIt;
   tip::TableRecord & output = *m_outputIt;
   tip::TableRecord & output2 = *m_output2It;

   dataSubselector::Cuts my_cuts(m_infile, m_evtable);
      
   my_cuts.addRangeCut("RA", "deg", 83, 93);
   my_cuts.addSkyConeCut(83., 22., 20);
   my_cuts.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                       dataSubselector::RangeCut::CLOSED, 1);
      
   tip::Index_t npts(0);
   tip::Index_t npts2(0);
      
   std::map<std::string, double> params;

   for (; m_inputIt != m_inputTable->end(); ++m_inputIt) {
      if (my_cuts.accept(input)) {
         output = input;
         ++m_outputIt;
         npts++;
      }
      input["ENERGY"].get(params["ENERGY"]);
      input["TIME"].get(params["TIME"]);
      input["RA"].get(params["RA"]);
      input["DEC"].get(params["DEC"]);
      params["CALIB_VERSION[1]"] = 1;
      if (my_cuts.accept(params)) {
         output2 = input;
         ++m_output2It;
         npts2++;
      }
      CPPUNIT_ASSERT(my_cuts.accept(input) == my_cuts.accept(params));

      if (my_cuts.accept(input)) {
         params["CALIB_VERSION[1]"] = 0;
         CPPUNIT_ASSERT(my_cuts.accept(input) != my_cuts.accept(params));
      }
   }

   m_outputTable->setNumRecords(npts);
   m_outputTable2->setNumRecords(npts2);
   
   my_cuts.writeDssKeywords(m_outputTable->getHeader());
   my_cuts.writeDssKeywords(m_outputTable2->getHeader());

   delete m_inputTable;
   delete m_outputTable;
   delete m_outputTable2;

   st_facilities::FitsUtil::writeChecksums(m_outfile);
   st_facilities::FitsUtil::writeChecksums(m_outfile2);

   dataSubselector::Cuts cuts1(m_outfile, m_evtable);
   dataSubselector::Cuts cuts2(m_outfile2, m_evtable);
   CPPUNIT_ASSERT(cuts1 == cuts2);

   dataSubselector::Cuts cuts(m_infile, m_evtable);
   CPPUNIT_ASSERT(!(cuts == cuts1));
}

void DssTests::compareCutsWithoutGtis() {
   dataSubselector::Cuts cuts1;
   cuts1.addRangeCut("Energy", "MeV", 30., 2e5);
   cuts1.addSkyConeCut(83.57, 22.01, 20);

   dataSubselector::Cuts cuts2(cuts1);

   dataSubselector::Gti gti1, gti2;
   gti1.insertInterval(100., 500.);
   gti2.insertInterval(300., 700.);

   cuts1.addGtiCut(gti1);
   cuts2.addGtiCut(gti2);

   CPPUNIT_ASSERT(!(cuts1 == cuts2));
   CPPUNIT_ASSERT(cuts1.compareWithoutGtis(cuts2));
   CPPUNIT_ASSERT(dataSubselector::Cuts::isTimeCut(cuts2[2]));
}

void DssTests::mergeGtis() {
   dataSubselector::Cuts cuts1;
   cuts1.addRangeCut("Energy", "MeV", 30., 2e5);
   cuts1.addSkyConeCut(83.57, 22.01, 20);

   dataSubselector::Cuts cuts2;
   
   cuts2 = cuts1;

   CPPUNIT_ASSERT(cuts2 == cuts1);

   dataSubselector::Gti gti1, gti2;
   gti1.insertInterval(100., 500.);
   gti2.insertInterval(300., 700.);

   cuts1.addGtiCut(gti1);
   cuts2.addGtiCut(gti2);

   std::vector<dataSubselector::Cuts> my_cuts;
   my_cuts.push_back(cuts1);
   my_cuts.push_back(cuts2);

   dataSubselector::Cuts newCuts = dataSubselector::Cuts::mergeGtis(my_cuts);

   CPPUNIT_ASSERT(newCuts != cuts1);
   CPPUNIT_ASSERT(newCuts != cuts2);
   CPPUNIT_ASSERT(newCuts.compareWithoutGtis(cuts1));
   CPPUNIT_ASSERT(newCuts.compareWithoutGtis(cuts2));

   std::vector<const dataSubselector::GtiCut *> gtiCuts;
   newCuts.getGtiCuts(gtiCuts);

   CPPUNIT_ASSERT(gtiCuts.size() == 1);

   dataSubselector::Gti mergedGti;
   mergedGti.insertInterval(100., 700.);

   CPPUNIT_ASSERT(!(gtiCuts.at(0)->gti() != mergedGti));
}

void DssTests::compareUnorderedCuts() {
   dataSubselector::Cuts cuts0;
   dataSubselector::Cuts cuts1;

   cuts0.addRangeCut("RA", "deg", 83, 93);
   cuts0.addSkyConeCut(83., 22., 20);
   cuts0.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                     dataSubselector::RangeCut::CLOSED, 1);
      
   cuts1.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                     dataSubselector::RangeCut::CLOSED, 1);
   cuts1.addSkyConeCut(83., 22., 20);
   cuts1.addRangeCut("RA", "deg", 83, 93);

   CPPUNIT_ASSERT(cuts0 == cuts1);
}

void DssTests::test_DssFormatting() {
   std::string testfile1("dss_test1.fits");
   std::string testfile2("dss_test2.fits");
   if (st_facilities::Util::fileExists(testfile1)) {
      std::remove(testfile1.c_str());
   }
   if (st_facilities::Util::fileExists(testfile2)) {
      std::remove(testfile2.c_str());
   }
   tip::IFileSvc::instance().createFile(testfile1, m_infile);
   tip::IFileSvc::instance().createFile(testfile2, m_infile);

   dataSubselector::Cuts cuts1;
//   cuts1.addRangeCut("TIME", "s", -1.5046090110e7, 505910.);
   cuts1.addRangeCut("ENERGY", "MeV", 30, 20000);
   tip::Table * table1 =
      tip::IFileSvc::instance().editTable(testfile1, "EVENTS");
   cuts1.writeDssKeywords(table1->getHeader());
   delete table1;

   tip::Table * table2 =
      tip::IFileSvc::instance().editTable(testfile2, "EVENTS");
   cuts1.writeDssKeywords(table2->getHeader());
   delete table2;

   table1 = tip::IFileSvc::instance().editTable(testfile1, "EVENTS");
   table2 = tip::IFileSvc::instance().editTable(testfile2, "EVENTS");

   const tip::Header & header1 = table1->getHeader();
   const tip::Header & header2 = table2->getHeader();

   std::string value1, value2;
   header1["DSVAL1"].get(value1);
   header2["DSVAL1"].get(value2);
   CPPUNIT_ASSERT(value1 == value2);

   delete table1;
   delete table2;

   CPPUNIT_ASSERT(!dataSubselector::Cuts::isTimeCut(cuts1[0]));
}

void DssTests::test_removeRangeCuts() {
   dataSubselector::Cuts cuts0;
   dataSubselector::Cuts cuts1;

   cuts0.addRangeCut("RA", "deg", 83, 93);
   cuts0.addSkyConeCut(83., 22., 20);
   cuts0.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                     dataSubselector::RangeCut::CLOSED, 1);
      
   cuts1.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                     dataSubselector::RangeCut::CLOSED, 1);
   cuts1.addSkyConeCut(83., 22., 20);


   CPPUNIT_ASSERT(cuts0 != cuts1);

   std::vector<dataSubselector::RangeCut *> removedCuts;
   cuts0.removeRangeCuts("RA", removedCuts);
   
   CPPUNIT_ASSERT(cuts0 == cuts1);

   dataSubselector::Cuts cuts2;
   for (size_t j = 0; j < removedCuts.size(); j++) {
      cuts2.addCut(*removedCuts.at(j));
   }
   
   dataSubselector::Cuts cuts3;
   cuts3.addRangeCut("RA", "deg", 83, 93);

   CPPUNIT_ASSERT(cuts2 == cuts3);
}

void DssTests::test_mergeRangeCuts() {
   dataSubselector::Cuts cuts0;

   cuts0.addRangeCut("ENERGY", "MeV", 100, 2e5);
   cuts0.addRangeCut("ENERGY", "MeV", 150, 2e5, 
                     dataSubselector::RangeCut::MINONLY);
   cuts0.addRangeCut("ENERGY", "MeV", 200, 1e5,
                     dataSubselector::RangeCut::MAXONLY);
   cuts0.addRangeCut("RA", "deg", 100, 200,
                     dataSubselector::RangeCut::MAXONLY, 1);
   cuts0.addRangeCut("RA", "deg", 0, 133,
                     dataSubselector::RangeCut::MAXONLY, 1);
   cuts0.addSkyConeCut(83., 22., 20);
   cuts0.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                      dataSubselector::RangeCut::CLOSED, 1);
   cuts0.mergeRangeCuts();

   dataSubselector::Cuts cuts1;
   cuts1.addRangeCut("ENERGY", "MeV", 150, 1e5);
   cuts1.addRangeCut("RA", "deg", 0, 133, 
                     dataSubselector::RangeCut::MAXONLY, 1);
   cuts1.addSkyConeCut(83., 22., 20);
   cuts1.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                      dataSubselector::RangeCut::CLOSED, 1);

   CPPUNIT_ASSERT(cuts0 == cuts1);
}

void DssTests::test_BitMaskCut() {
   dataSubselector::BitMaskCut cut("EVENT_CLASS", 4);
   std::map<std::string, double> pars;
   pars["EVENT_CLASS"] = 4;
   CPPUNIT_ASSERT(cut.accept(pars));
   pars["EVENT_CLASS"] = 5;
   CPPUNIT_ASSERT(cut.accept(pars));
   pars["EVENT_CLASS"] = 6;
   CPPUNIT_ASSERT(cut.accept(pars));
   pars["EVENT_CLASS"] = 7;
   CPPUNIT_ASSERT(cut.accept(pars));
   pars["EVENT_CLASS"] = 8;
   CPPUNIT_ASSERT(!cut.accept(pars));
   CPPUNIT_ASSERT(cut.filterString() == "((EVENT_CLASS/4)%2 == 1)");

   /// Check post-Pass 7 behavior.
   // 0324 (octal) = 212 (decimal) = 11010100 (binary)
   dataSubselector::BitMaskCut p8_cut("EVENT_TYPE", 0324, "P8");
   std::string filter(p8_cut.filterString());
   CPPUNIT_ASSERT(filter == "((EVENT_TYPE&o324) != o0)");

   pars["EVENT_TYPE"] = 212;
   CPPUNIT_ASSERT(p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 212 + 1;
   CPPUNIT_ASSERT(p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 212 + 2;
   CPPUNIT_ASSERT(p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 212 + 8;
   CPPUNIT_ASSERT(p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 212 + 32;
   CPPUNIT_ASSERT(p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 1;
   CPPUNIT_ASSERT(!p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 2;
   CPPUNIT_ASSERT(!p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 3;
   CPPUNIT_ASSERT(!p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 8;
   CPPUNIT_ASSERT(!p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 10;
   CPPUNIT_ASSERT(!p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 11;
   CPPUNIT_ASSERT(!p8_cut.accept(pars));

   pars["EVENT_TYPE"] = 256;
   CPPUNIT_ASSERT(!p8_cut.accept(pars));

   /// Check if a candidate cut supercedes the current cut using the
   /// validity mask mechanism.
   std::string evclassFile("evclass_validity_masks.txt");
   std::ofstream outfile1(evclassFile.c_str());
   outfile1 << "128,1920\n"
            << "256,1792\n"
            << "512,1536\n";
   outfile1.close();

   std::string evtypeFile("evtype_validity_masks.txt");
   std::ofstream outfile2(evtypeFile.c_str());
   outfile2 << "3,1023\n"
            << "60,1023\n"
            << "960,1023\n"
            << "56,56\n"
            << "1,1\n"
            << "832,832\n";
   outfile2.close();

   dataSubselector::BitMaskCut::setValidityMasks(evclassFile, evtypeFile);
   std::remove(evclassFile.c_str());
   std::remove(evtypeFile.c_str());

   dataSubselector::BitMaskCut source("EVENT_CLASS", 128, "P8");
   dataSubselector::BitMaskCut clean("EVENT_CLASS", 256, "P8");
   dataSubselector::BitMaskCut ultraclean("EVENT_CLASS", 512, "P8");

   CPPUNIT_ASSERT(ultraclean.supercedes(source));
   CPPUNIT_ASSERT(ultraclean.supercedes(clean));
   CPPUNIT_ASSERT(clean.supercedes(source));

   try {
      source.supercedes(clean);
   } catch (std::runtime_error &) {
   }
   try {
      source.supercedes(ultraclean);
   } catch (std::runtime_error &) {
   } 

   dataSubselector::BitMaskCut fb("EVENT_TYPE", 3, "P8");
   dataSubselector::BitMaskCut front("EVENT_TYPE", 1, "P8");
   dataSubselector::BitMaskCut psf("EVENT_TYPE", 60, "P8");
   dataSubselector::BitMaskCut edisp("EVENT_TYPE", 960, "P8");
   dataSubselector::BitMaskCut psf123("EVENT_TYPE", 56, "P8");
   dataSubselector::BitMaskCut edisp023("EVENT_TYPE", 832, "P8");
   
   CPPUNIT_ASSERT(psf.supercedes(fb));
   CPPUNIT_ASSERT(fb.supercedes(psf));
   CPPUNIT_ASSERT(edisp.supercedes(psf));

   CPPUNIT_ASSERT(psf123.supercedes(fb));
   CPPUNIT_ASSERT(psf123.supercedes(psf));
   CPPUNIT_ASSERT(psf123.supercedes(edisp));

   CPPUNIT_ASSERT(edisp023.supercedes(fb));
   CPPUNIT_ASSERT(edisp023.supercedes(psf));
   CPPUNIT_ASSERT(edisp023.supercedes(edisp));

   try {
      psf123.supercedes(edisp023);
   } catch (std::runtime_error &) {
   }
   try {
      edisp023.supercedes(psf123);
   } catch (std::runtime_error &) {
   }
   try {
      psf.supercedes(psf123);
   } catch (std::runtime_error &) {
   }
   try {
      edisp.supercedes(edisp023);
   } catch (std::runtime_error &) {
   }
   try {
      fb.supercedes(psf123);
   } catch (std::runtime_error &) {
   }
   try {
      fb.supercedes(edisp023);
   } catch (std::runtime_error &) {
   }

   dataSubselector::BitMaskCut invalid("EVENT_TYPE", 1024, "P8");
   try {
      psf.supercedes(invalid);
      CPPUNIT_ASSERT(false);  // runtime_error exception should be thrown.
   } catch (std::runtime_error &) {
      CPPUNIT_ASSERT(true);
   }
}

void DssTests::test_VersionCut() {
   dataSubselector::VersionCut cut("IRF_VERSION", "V6MC");
   std::map<std::string, double> pars;  
   // Can leave pars empty since all accept methods should return true.
   CPPUNIT_ASSERT(cut.accept(pars));

   CPPUNIT_ASSERT(cut.colname() == "IRF_VERSION");
   CPPUNIT_ASSERT(cut.version() == "V6MC");
   CPPUNIT_ASSERT(cut.version() != "V7");

   dataSubselector::VersionCut new_cut("IRF_VERSION", "V7");
   CPPUNIT_ASSERT(!(new_cut == cut));
   CPPUNIT_ASSERT(new_cut.supercedes(cut));
   CPPUNIT_ASSERT(new_cut.type() == "version");

   dataSubselector::VersionCut cut_copy(new_cut);
   CPPUNIT_ASSERT(cut_copy.type() == "version");
}

void DssTests::test_irfName() {
   dataSubselector::Cuts cuts1;
   cuts1.setIrfs("P7SOURCE_V6::BACK");
   CPPUNIT_ASSERT(cuts1.bitMaskCut()->mask() == 4);
   try {
      cuts1.setIrfs("P7TRANSIENT_V6");
   } catch (std::runtime_error & eObj) {
      st_facilities::Util::expectedException(eObj, "Bit mask cut already set");
   }

   CPPUNIT_ASSERT(cuts1.conversionTypeCut()->minVal() == 1);
   CPPUNIT_ASSERT(cuts1.conversionTypeCut()->maxVal() == 1);
   CPPUNIT_ASSERT(cuts1.pass_ver() == "P7V6");

   dataSubselector::Cuts cuts2;
   cuts2.setIrfs("P7TRANSIENT_V6");
   CPPUNIT_ASSERT(cuts2.bitMaskCut()->mask() == 1);
   CPPUNIT_ASSERT(cuts2.conversionTypeCut() == 0);
   CPPUNIT_ASSERT(cuts2.pass_ver() == "P7V6");

   dataSubselector::Cuts cuts3;
   cuts3.setIrfs("P7REP_SOURCE_V10");
   CPPUNIT_ASSERT(cuts3.bitMaskCut()->mask() == 4);
   CPPUNIT_ASSERT(cuts3.conversionTypeCut() == 0);
   CPPUNIT_ASSERT(cuts3.pass_ver() == "P7REP");

   dataSubselector::Cuts cuts4;
   cuts4.setIrfs("P7SOURCE_V6");
   CPPUNIT_ASSERT(cuts4.CALDB_implied_irfs() == "P7SOURCE_V6");

   dataSubselector::Cuts cuts5;
   cuts5.setIrfs("P7REP_SOURCE_V10");
   CPPUNIT_ASSERT(cuts5.CALDB_implied_irfs() != "P7REP_SOURCE_V10");
}

void DssTests::test_rangeCut() {
   dataSubselector::Cuts my_cuts;
   my_cuts.addRangeCut("CALIB_VERSION", "dimensionless", 1, 1,
                       dataSubselector::RangeCut::CLOSED);
   std::map<std::string, double> params;
   params["CALIB_VERSION"] = 1;
   CPPUNIT_ASSERT(my_cuts.accept(params));
}

int main(int iargc, char * argv[]) {

   if (iargc > 1 && std::string(argv[1]) == "-d") {
      DssTests testObj;
      testObj.setUp();
      testObj.test_BitMaskCut();
      testObj.tearDown();
   } else {
      CppUnit::TextTestRunner runner;
      runner.addTest(DssTests::suite());
      bool result = runner.run();
      if (result) {
         return 0;
      } else {
         return 1;
      }
   }
}
