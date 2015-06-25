/**
 * @file test.cxx
 * @brief Test program for st_facilities
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/st_facilities/src/test/test.cxx,v 1.14 2009/06/26 21:34:52 jchiang Exp $
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "st_facilities/dgaus8.h"
#include "st_facilities/GaussianQuadrature.h"
#include "PowerLaw.h"

#include "st_facilities/Env.h"
#include "st_facilities/FileSys.h"
#include "st_facilities/FitsImage.h"
#include "st_facilities/Util.h"

using namespace st_facilities;

class st_facilitiesTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(st_facilitiesTests);

   CPPUNIT_TEST(test_dgaus8);
   CPPUNIT_TEST(test_GaussianQuadrature);
   CPPUNIT_TEST_EXCEPTION(test_Util_file_ok, std::runtime_error);
   CPPUNIT_TEST(test_Util_readLines);
   CPPUNIT_TEST(test_Util_expectedException);
   CPPUNIT_TEST(test_Util_resolve_fits_files);
   CPPUNIT_TEST(test_Env_appendNames);
//   CPPUNIT_TEST(test_Env_expandEnvVar);
//    CPPUNIT_TEST(test_Env_getDataDir);
//    CPPUNIT_TEST(test_FileSys_expandFileList);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();

   void test_dgaus8();
   void test_GaussianQuadrature();
   void test_Util_file_ok();
   void test_Util_readLines();
   void test_Util_expectedException();
   void test_Util_resolve_fits_files();
   void test_Env_appendNames();
   void test_Env_expandEnvVar();
   void test_Env_getDataDir();
   void test_FileSys_expandFileList();

private:

   std::string m_filename;

};

void st_facilitiesTests::setUp() {
   m_filename = "test_file.txt";
   std::ofstream file(m_filename.c_str());
   file << "line 0\n"
        << "line 1\n"
        << "line 2\n"
        << "#line 3\n"
        << " "
        << "\n";
   file.close();
}

void st_facilitiesTests::tearDown() {
   std::remove(m_filename.c_str());
}

class Linear {
public:
   double operator()(double x) const {
      return x;
   }
};

class Foo {
public:
   Foo(const Linear & bar) : m_bar(bar) {}
   double operator()(double x) const {
      double err;
      int ier;
      return x*GaussianQuadrature::dgaus8(m_bar, 0, 1, err, ier);
   }
private:
   const Linear & m_bar;
};

class Edisp {
public:
   Edisp(double ltail, double rwidth) : m_ltail(ltail), m_rwidth(rwidth) {}
   double operator()(double x) const {
      double arg(x/m_rwidth);
      if (arg > 40) {
         return std::pow(x + 1, m_ltail)*std::exp(-arg);
      }
      return std::pow(x + 1, m_ltail)/(1 + std::exp(arg));
   }
private:
   double m_ltail;
   double m_rwidth;
};

void st_facilitiesTests::test_dgaus8() {
   double err(1e-5);
   double result(0);
   int ier;
   double tol(1e-4);

   PowerLaw pl(1, -2);
   result = GaussianQuadrature::dgaus8(pl, 1, 5, err, ier);
   double trueValue(0.8);
   CPPUNIT_ASSERT(std::fabs((result - trueValue)/trueValue) < tol);

   Linear fx;
   result = GaussianQuadrature::dgaus8(fx, 0, 1, err, ier);
   trueValue = 0.5;
   CPPUNIT_ASSERT(std::fabs((result - trueValue)/trueValue) < tol);

   Foo fxy(fx);
   result = GaussianQuadrature::dgaus8(fxy, 0, 1, err, ier);
   trueValue = 0.25;
   CPPUNIT_ASSERT(std::fabs((result - trueValue)/trueValue) < tol);

   double ltail(20);
   double rwidth(0.1);

   double xmin(-1);
   double xmax(2);

   Edisp edisp(ltail, rwidth);

   result = GaussianQuadrature::dgaus8(edisp, xmin, xmax, err, ier);

   size_t npts(500);
   std::vector<double> xx;
   std::vector<double> yy;

   double dx((xmax - xmin)/(npts-1));
   for (size_t i(0); i < npts; i++) {
      xx.push_back(i*dx + xmin);
      yy.push_back(edisp(xx.back()));
   }
   
   double integral(0);
   for (size_t i(1); i < npts; i++) {
      integral += (yy.at(i) + yy.at(i-1))/2.*(xx.at(i) - xx.at(i-1));
   }

   CPPUNIT_ASSERT(std::fabs((result - integral)/integral) < tol);
}

PowerLaw powerLaw(1., 2.);

double power_law(double * x) {
   return powerLaw(*x);
}

void st_facilitiesTests::test_GaussianQuadrature() {
   double xmin(0);
   double xmax(4.);
   double err(1e-5);
   long ier;

   double true_value;
   if (powerLaw.index() != 1.) {
      true_value = powerLaw.prefactor()*(pow(xmax, powerLaw.index() + 1.)
                                         - pow(xmin, powerLaw.index() + 1.))
                                         /(powerLaw.index() + 1.);
   } else {
      true_value = powerLaw.prefactor()*log(xmax/xmin);
   }
   double tol(1e-4);

   double result(GaussianQuadrature::integrate(&power_law, xmin, xmax,
                                               err, ier));

   CPPUNIT_ASSERT(std::fabs((result - true_value)/true_value) < tol);
}

void st_facilitiesTests::test_Util_file_ok() {
   std::string filename("foo");
   std::remove(filename.c_str());
   Util::file_ok(filename);
}

void st_facilitiesTests::test_Util_readLines() {
   std::vector<std::string> lines;
   Util::readLines(m_filename, lines);
   CPPUNIT_ASSERT(lines.size() == 3);
   Util::readLines(m_filename, lines, "");
   CPPUNIT_ASSERT(lines.size() == 4);
}

void st_facilitiesTests::test_Util_expectedException() {
   try {
      test_Util_file_ok();
   } catch (std::exception & eObj) {
      CPPUNIT_ASSERT(Util::expectedException(eObj, "File not found"));
      CPPUNIT_ASSERT(!Util::expectedException(eObj, "File not fund"));
   }
}

void st_facilitiesTests::test_Util_resolve_fits_files() {
   std::vector<std::string> lines;
   Util::resolve_fits_files(m_filename, lines);
   CPPUNIT_ASSERT(lines.size() == 3);

   Util::resolve_fits_files("@" + m_filename, lines);
   CPPUNIT_ASSERT(lines.size() == 3);
}

void st_facilitiesTests::test_Env_appendNames() {
   std::string tmp_str;
   std::string::size_type pos;
   std::string seg1 = "seg1-name";
   std::string seg2 = "seg2-name";

   tmp_str = Env::appendPath("", "");
   CPPUNIT_ASSERT(tmp_str.empty());

   tmp_str = Env::appendPath(seg1, "");
   CPPUNIT_ASSERT(tmp_str == seg1);

   tmp_str = Env::appendPath("", seg2);
   CPPUNIT_ASSERT(tmp_str == seg2);

   tmp_str = Env::appendPath(seg1, seg2);
   pos = tmp_str.find(seg1);
   CPPUNIT_ASSERT(0 == pos);
   pos = tmp_str.find(seg2);
   CPPUNIT_ASSERT(tmp_str.size() - seg2.size() == pos);
}

void st_facilitiesTests::test_Env_expandEnvVar() {
   std::string to_expand;
   std::string expanded;

   to_expand = "$STTEST";
   Env::expandEnvVar(to_expand, expanded);
   CPPUNIT_ASSERT("sttest" == expanded);
   expanded.clear();

   to_expand = "${STTEST}";
   Env::expandEnvVar(to_expand, expanded);
   CPPUNIT_ASSERT(expanded == "sttest");
   expanded.clear();

   to_expand = "$(STTEST)";
   Env::expandEnvVar(to_expand, expanded);
   CPPUNIT_ASSERT(expanded == "sttest");
   expanded.clear();

   to_expand = "$(STTEST) filler $STTEST";
   Env::expandEnvVar(to_expand, expanded);
   CPPUNIT_ASSERT(expanded == "sttest filler sttest");
   expanded.clear();

   to_expand = "filler$(STTEST)filler$STTEST";
   Env::expandEnvVar(to_expand, expanded);
   CPPUNIT_ASSERT(expanded == "fillersttestfillersttest");
   expanded.clear();

   to_expand = "$(STTEST)}";
   Env::expandEnvVar(to_expand, expanded);
   CPPUNIT_ASSERT(expanded == "sttest}");
   expanded.clear();

   to_expand = "$STTEST filler";
   Env::expandEnvVar(to_expand, expanded);
   CPPUNIT_ASSERT(expanded == "sttest filler");
   expanded.clear();

   to_expand = "$";
   try {
      Env::expandEnvVar(to_expand, expanded);
      CPPUNIT_ASSERT(false);
   } catch (const std::exception &) {
      CPPUNIT_ASSERT(expanded == to_expand);
   }
   expanded.clear();

   to_expand = "$$";
   try {
      Env::expandEnvVar(to_expand, expanded);
      CPPUNIT_ASSERT(false);
   } catch (const std::exception &) {
      CPPUNIT_ASSERT(expanded == to_expand);
   }
   expanded.clear();

   to_expand = "$STTEST_run_on";
   try {
      Env::expandEnvVar(to_expand, expanded);
      CPPUNIT_ASSERT(false);
   } catch (const std::exception &) {
      CPPUNIT_ASSERT(expanded == to_expand);
   }
   expanded.clear();

   to_expand = "${STTEST:f}";
   try {
      Env::expandEnvVar(to_expand, expanded);
      CPPUNIT_ASSERT(false);
   } catch (const std::exception &) {
      CPPUNIT_ASSERT(expanded == to_expand);
   }
   expanded.clear();

   to_expand = "prefix${STTEST";
   try {
      Env::expandEnvVar(to_expand, expanded);
      CPPUNIT_ASSERT(false);
   } catch (const std::exception &) {
      CPPUNIT_ASSERT(expanded == to_expand);
   }
   expanded.clear();
}

void st_facilitiesTests::test_Env_getDataDir() {
   // For comparison of output, need local and install areas.
   std::string install_area;
   std::string local_area = Env::getEnv("ST_FACILITIESROOT");

   try {
      // For comparison purposes, expand the install area env variable.
      std::string install_area = Env::getEnv("DATADIR");
   } catch (const std::exception &) {
      // That's OK.
   }

   std::string data_dir;
   // First test: an invalid package identifier.
   try {
      data_dir = Env::getDataDir("invalid_pkg");
      // If the above line did not throw an exception, the data_dir
      // must be the install area.
      CPPUNIT_ASSERT(!install_area.empty() && install_area == data_dir);
   } catch (const std::exception &) {
      // Expected only if install area env variable is not set.
      CPPUNIT_ASSERT(install_area.empty());
   }

   data_dir.erase();
   try {
      data_dir = Env::getDataDir("st_facilities");
      // If the above line did not throw an exception, the data_dir
      // must be the local area.
      CPPUNIT_ASSERT(0 == data_dir.find(local_area));
   } catch (const std::exception &) {
      // Unexpected.
      CPPUNIT_ASSERT(false);
   }
}

void st_facilitiesTests::test_FileSys_expandFileList() {
  std::string list_file;
  FileSys::FileNameCont cont;

  // Expansion of a file which doesn't exist should throw.
  try {
    list_file = "@a-non-existent-file";
    cont = FileSys::expandFileList(list_file);
    CPPUNIT_ASSERT(false);
  } catch (const std::exception &) {
  }

  // Non-expansion case: no leading @.
  list_file = Env::appendFileName(Env::appendFileName("$ST_FACILITIESROOT",
                                                      "data"), "list_file");
  cont = FileSys::expandFileList(list_file);
  // Container should hold just the original file name because no
  // expansion occurred.
  CPPUNIT_ASSERT(1 == cont.size());
  std::string expanded_list;
  Env::expandEnvVar(list_file, expanded_list);
  CPPUNIT_ASSERT(*cont.begin() == expanded_list);

  // Expansion of a file which exists and contains a list of files
  // should be done correctly.  The test file contains itself, as well
  // as a second line with an arbitrary name.
  list_file = "@" + list_file;
  cont = FileSys::expandFileList(list_file);
  CPPUNIT_ASSERT(2 == cont.size());

// The first line in list_file contains non-Windows path separators; therefore
// the following assert will fail on Windows.
//  CPPUNIT_ASSERT(*cont.begin() == expanded_list);
  CPPUNIT_ASSERT(cont.back() == "fits_file1.fits");
}

int main() {
   CppUnit::TextTestRunner runner;

   runner.addTest(st_facilitiesTests::suite());

   bool result = runner.run();
   if (result) {
      return 0;
   } else {
      return 1;
   }
}
