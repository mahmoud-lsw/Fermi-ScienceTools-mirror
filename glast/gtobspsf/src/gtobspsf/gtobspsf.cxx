// $Id: gtobspsf.cxx,v 1.11 2009/10/07 18:00:15 elwinter Exp $

// N.B. I am currently using abort() to handle run-time exceptions,
// since I don't thoroughly understand how to use C++ exceptions yet.

// N.B. I am making a lot of errors with container methods, so I
// dropped back to using calloc() for now.

// This file contains the main code for gtobspsf, which computes an
// observed Point Spread Function (PSF) from a GLAST event (FT1)
// file. The output from this program is intended to be compatible
// with the output from gtpsf, in that the angular integral of any
// single PSF should sum to approximately unity.

// The computation is performed in the following steps:

// 1. Bin each event by energy; energy bins are equally-spaced in logE
//    space.

// 2. Within each energy bin, sub-bin each event by angle relative to
//    the source location.

// 3. Compute the error for each bin using the "Gehrels fudge".

// 4. Read the total exposure time from the ONTIME keyword of the GTI
//    table. Use it to compute the total exposure as effarea*ontime.

// 5. Compute the uncorrected intensity in each bin by dividing the
//    counts (and errors) in each angular bin by the angular area of the
//    bin, and the total exposure. These values are in units of
//    photons/sr.

// 6. Subtract the constant background intensity from the intensity for
//    each bin. Any negative intensities are left unchanged.

// 7. Compute the PSF in units of sr^-1 by multiplying the corrected
//    intensities and errors by the total exposure, and dividing by
//    the total number of photons in all bins.

//*****************************************************************************

// Standard headers
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <memory>

// Other SAE headers
#include "astro/SkyDir.h"
#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"
#include "st_facilities/Util.h"
#include "st_stream/StreamFormatter.h"
#include "tip/IFileSvc.h"
#include "tip/Header.h"
#include "tip/Table.h"

// Project headers

//*****************************************************************************

// Local macros

// Conversion factors between radians and degrees.
#define DEGREES_PER_RADIAN (180.0 / M_PI)
#define RADIANS_PER_DEGREE (M_PI / 180.0)

// 2 pi
#define TWO_PI (2 * M_PI)

//*****************************************************************************

class GtobspsfApp : public st_app::StApp {

public:

  // Define the required run() method.
  virtual void run();

  // ?? This seems to be required at compile time.
  virtual ~GtobspsfApp() throw() {}

private:

  // Debugging output flag.
  bool m_debug;

  // Number of energy bins.
  int m_nenergies;

  // Array of lower-bound energies (MeV) for logarithmic energy bins.
  double *m_bin_energies;

  // Number of angle bins.
  int m_nangles;

  // Array of lower-bound angles (degrees) for angle bins.
  double *m_bin_angles;

  // Array of (m_nenergies x m_nangles) values for the normalized PSF
  // at each energy.
  double **m_PSF;

  // Array of (m_nenergies x m_nangles) values for the normalized PSF
  // errors at each energy.
  double **m_PSF_error;

  // Constant exposure (cm^2 s) for all bins.
  double m_exposure;

  // Method to create output FITS file.
  void writeFitsFile();

};

//*****************************************************************************

void GtobspsfApp::run() {

  //---------------------------------------------------------------------------

  // Fetch and validate run parameters.

  // Get the group of parameters for this application.
  st_app::AppParGroup & pars(getParGroup("gtobspsf"));

  // Check if debugging mode was requested.
//   bool debug;
//   if (pars["debug"] == "yes") {
//     debug = true;
//   } else {
//     debug = false;
//   }
//   m_debug = debug;
  m_debug = false;

  // Event (FT1) file name
  pars.Prompt("evfile");
  std::string evfile = pars["evfile"];

  // Event table name
  std::string evtable = pars["evtable"];

  // Constant detector effective area
  pars.Prompt("effarea");
  double effarea = pars["effarea"];
  if (effarea <= 0.0) {
    std::cout << "effarea parameter must be > 0!" << std::endl;
    abort();
  }

  // Constant background flux
  pars.Prompt("bgflux");
  double bgflux = pars["bgflux"];
  if (bgflux < 0.0) {
    std::cout << "bgflux parameter must be >= 0!" << std::endl;
    abort();
  }

  // Output file name for PSF results
  pars.Prompt("outfile");
  std::string outfile = pars["outfile"];

  // Output PSF extension name
  std::string outtable = pars["outtable"];

  // Source right ascension (degrees)
  pars.Prompt("src_ra");
  double src_ra = pars["src_ra"];
  if (src_ra < 0.0 || src_ra > 360.0) {
    std::cout << "src_ra must be >= 0 and <= 360!" << std::endl;
    abort();
  }

  // Source declination (degrees)
  pars.Prompt("src_dec");
  double src_dec = pars["src_dec"];
  if (src_dec < -90.0 || src_dec > 90.0) {
    std::cout << "src_dec must be >= -90 and <= 90!" << std::endl;
    abort();
  }

  // Maximum off-axis angle (degrees) to consider for PSF
  pars.Prompt("maxangle");
  double maxangle = pars["maxangle"];
  if (maxangle <= 0.0 || maxangle > 180.0) {
    std::cout << "maxangle must be > 0 and <= 180!" << std::endl;
    abort();
  }

  // Number of equal-size angular bins for PSF
  pars.Prompt("nangles");
  int nangles = pars["nangles"];
  if (nangles <= 0) {
    std::cout << "nangles must be > 0!" << std::endl;
    abort();
  }
  m_nangles = nangles;

  // Minimum energy value for log-spaced energy bins (MeV).
  pars.Prompt("emin");
  double emin = pars["emin"];
  if (emin <= 0.0) {
    std::cout << "emin must be > 0!" << std::endl;
    abort();
  }

  // Maximum energy value for log-spaced energy bins (MeV).
  pars.Prompt("emax");
  double emax = pars["emax"];
  if (emax <= 0.0 || emax <= emin) {
    std::cout << "emax must be > emin!" << std::endl;
    abort();
  }

  // Number of log-spaced energy bins.
  pars.Prompt("nenergies");
  int nenergies = pars["nenergies"];
  if (nenergies <= 0) {
    std::cout << "nenergies must be > 0!" << std::endl;
    abort();
  }
  m_nenergies = nenergies;

  // Save the current parameter values.
  pars.Save();

  if (m_debug) {
    std::cout << "evfile = " << evfile << std::endl;
    std::cout << "evtable = " << evtable << std::endl;
    std::cout << "effarea = " << effarea << std::endl;
    std::cout << "bgflux = " << bgflux << std::endl;
    std::cout << "outfile = " << outfile << std::endl;
    std::cout << "outtable = " << outtable << std::endl;
    std::cout << "src_ra = " << src_ra << std::endl;
    std::cout << "src_dec = " << src_dec << std::endl;
    std::cout << "maxangle = " << maxangle << std::endl;
    std::cout << "nangles = " << nangles << std::endl;
    std::cout << "emin = " << emin << std::endl;
    std::cout << "emax = " << emax << std::endl;
    std::cout << "nenergies = " << nenergies << std::endl;
  }

  //---------------------------------------------------------------------------

  // Compute the (log) size of the energy bins.
  double dlogE = log(emax / emin) / m_nenergies;
  if (m_debug) {
    std::cout << "dlogE = " << dlogE << std::endl;
  }

  // Compute the lower bounds for each log-spaced energy bin.
  m_bin_energies = (double *) calloc(m_nenergies, sizeof(double));
  for (int i = 0; i < m_nenergies; ++i) {
    m_bin_energies[i] = emin * exp(i * dlogE);
  }
  if (m_debug) {
    for (int i = 0; i < m_nenergies; ++i) {
      std::cout << "m_bin_energies[" << i << "] = " << m_bin_energies[i] <<
	std::endl;
    }
  }

  //---------------------------------------------------------------------------

  // Compute the size of the angle bins.
  double dtheta = maxangle / m_nangles;
  if (m_debug) {
    std::cout << "dtheta = " << dtheta << std::endl;
  }

  // Compute the lower limit of each angle bin.
  m_bin_angles = (double *) calloc(m_nangles, sizeof(double));
  for (int j = 0; j < m_nangles; ++j) {
    m_bin_angles[j] = j * dtheta;
  }
  if (m_debug) {
    for (int j = 0; j < m_nangles; ++j) {
      std::cout << "m_bin_angles[" << j << "] = " << m_bin_angles[j] <<
	std::endl;
    }
  }

  //---------------------------------------------------------------------------

  // Bin the events by energy and angle.

  // Create an (m_nenergies x m_nangles) array to hold the counts in
  // each energy and angle bin.
  int **bins = (int **) calloc(m_nenergies, sizeof(int *));
  for (int i = 0; i < m_nenergies; ++i) {
    bins[i] = (int *) calloc(m_nangles, sizeof(int));
  }

  // Create a sky direction object for the center location for the PSF
  // computation (arguments are in degrees).
  astro::SkyDir center_direction(src_ra, src_dec);

  // Open the event table for read access.
  std::auto_ptr<const tip::Table>
    event(tip::IFileSvc::instance().readTable(evfile, evtable));

  // Bin each event by energy and angle.
  int n_events = 0;
  int n_good_events = 0;
  for (tip::Table::ConstIterator record = event->begin();
       record != event->end(); ++record) {

    // Increment the event count.
    ++n_events;
    if (m_debug) {
      std::cout << "n_events = " << n_events << std::endl;
    }

    // For clarity, make a local reference to the record contained in
    // the record iterator.
    tip::Table::ConstRecord & event(*record);

    // Get the position and energy fields from the event record.
    double ra = event["RA"].get();
    double dec = event["DEC"].get();
    double energy = event["ENERGY"].get();
    if (m_debug) {
      std::cout << "ra = " << ra << ", dec = " << dec << ", energy = " <<
	energy << std::endl;
    }

    // Skip this event if it is outside of the specified energy range.
    if (energy < emin || energy >= emax) {
      if (m_debug) {
	std::cout << "Skipping event " << n_events <<
	  " since energy is out of range." << std::endl;
      }
      continue;
    }

    // Create a sky direction object for the event position.
    astro::SkyDir direction(ra, dec);

    // Compute the angle (in degrees) from the source position to the
    // event position.
    double theta = center_direction.difference(direction) * DEGREES_PER_RADIAN;
    if (m_debug) {
      std::cout << "theta = " << theta << std::endl;
    }

    // Skip this event if it is outside of the specified angle range.
    if (theta >= maxangle) {
      if (m_debug) {
	std::cout << "Skipping event " << n_events <<
	  " since off-axis angle is too high." << std::endl;
      }
      continue;
    }

    // Increment the count of good events.
    ++n_good_events;
    if (m_debug) {
      std::cout << "n_good_events = " << n_good_events << std::endl;
    }

    // Compute the index of the log(E) bin.
    int i_energy = (int) floor(log(energy / emin) / dlogE);
    if (m_debug) {
      std::cout << "i_energy = " << i_energy << std::endl;
    }

    // Compute the index of the angle bin.
    int j_angle = (int) floor(theta / dtheta);

    // Increment the event count for this bin.
    bins[i_energy][j_angle]++;

  }
  if (m_debug) {
    for (int i = 0; i < m_nenergies; ++i) {
      for (int j = 0; j < m_nangles; ++j) {
 	std::cout << "bins[" << i << "][" << j << "] = " << bins[i][j] <<
	  std::endl;
      }
    }
    std::cout << "n_good_events = " << n_good_events << std::endl;
  }

  //---------------------------------------------------------------------------

  // Compute the error in each bin count using the "Gehrels fudge":
  // err = 1 + sqrt(n + 0.75)

  // Create the (m_nenergies x m_nangles) array to hold the error
  // extimates for each bin.
  double **bin_errors = (double **) calloc(m_nenergies, sizeof(double *));
  for (int i = 0; i < m_nenergies; ++i) {
    bin_errors[i] = (double *) calloc(m_nangles, sizeof(double));
  }

  // Estimate the error in each bin.
  for (int i = 0; i < m_nenergies; ++i) {
    for (int j = 0; j < m_nangles; ++j) {
      bin_errors[i][j] = 1.0 + sqrt(bins[i][j] + 0.75);
    }
  }
  if (m_debug) {
    for (int i = 0; i < m_nenergies; ++i) {
      for (int j = 0; j < m_nangles; ++j) {
 	std::cout << "bin_errors[" << i << "][" << j << "] = " <<
	  bin_errors[i][j] << std::endl;
      }
    }
  }

  //---------------------------------------------------------------------------

  // Fetch the total exposure time from the GTI table.

  // Open the GTI table for read access.
  std::auto_ptr<const tip::Table>
    gti(tip::IFileSvc::instance().readTable(evfile, "GTI"));

  // Get the GTI table header.
  const tip::Header & header(gti->getHeader());

  // Fetch the value of the ONTIME keyword.
  double ontime;
  header["ONTIME"].get(ontime);
  if (m_debug) {
    std::cout << "ontime = " << std::setprecision(12) << ontime << std::endl;
  }

  // Compute the total exposure.
  m_exposure = effarea * ontime;
  if (m_debug) {
    std::cout << "m_exposure = " << m_exposure << std::endl;
  }

  //---------------------------------------------------------------------------

  // Compute the PSF (sr^-1).

  // Compute the angular area (in steradians) of each angular
  // bin. Each bin is an annulus around the target point, of angular
  // width dtheta (degrees).
  double *domega = (double *) calloc(m_nangles, sizeof(double));
  for (int j = 0; j < m_nangles; ++j) {
    domega[j] = 2.0 * M_PI *
      (cos(j * dtheta * RADIANS_PER_DEGREE) -
       cos((j + 1) * dtheta * RADIANS_PER_DEGREE));
  }
  if (m_debug) {
    for (int j = 0; j < m_nangles; ++j) {
      std::cout << "domega[" << j << "] = " << domega[j] << std::endl;
    }
  }

  // Convert the bin counts and errors to intensities by dividing by
  // the angular area of the bin, the effective area of the detector,
  // and the total exposure time.
  double **intensity = (double **) calloc(m_nenergies, sizeof(double *));
  double **intensity_error = (double **) calloc(m_nenergies, sizeof(double *));
  for (int i = 0; i < m_nenergies; ++i) {
    intensity[i] = (double *) calloc(m_nangles, sizeof(double));
    intensity_error[i] = (double *) calloc(m_nangles, sizeof(double));
  }
  for (int i = 0; i < m_nenergies; ++i) {
    for (int j = 0; j < m_nangles; ++j) {
      intensity[i][j] = (double) bins[i][j] / domega[j] / m_exposure;
      intensity_error[i][j] = bin_errors[i][j] / domega[j] / m_exposure;
    }
  }
  if (m_debug) {
    for (int i = 0; i < m_nenergies; ++i) {
      for (int j = 0; j < m_nangles; ++j) {
   	std::cout << "intensity[" << i << "][" << j << "] = " <<
 	  intensity[i][j] << std::endl;
     	std::cout << "intensity_error[" << i << "][" << j  << "] = " <<
     	  intensity_error[i][j] << std::endl;
      }
    }
  }

  // N.B. At this point, all values of intensity and intensity_error
  // must be >= 0.

  // Subtract the constant background intensity (photons/cm^2/s/sr).
  double background_intensity = bgflux / TWO_PI;
  for (int i = 0; i < m_nenergies; ++i) {
    for (int j = 0; j < m_nangles; ++j) {
      intensity[i][j] -= background_intensity;
    }
  }
  if (m_debug) {
    for (int i = 0; i < m_nenergies; ++i) {
      for (int j = 0; j < m_nangles; ++j) {
  	std::cout << "intensity[" << i << "][" << j << "] = " <<
 	  intensity[i][j] << std::endl;
      }
    }
  }

  // N.B. Now that background has been subtracted, the intensities may
  // be negative. Leave them that way.

  // Convert the background-corrected intensities back to units of
  // sr^-1. This is the form that the PSF should be saved in.
  m_PSF = (double **) calloc(m_nenergies, sizeof(double *));
  m_PSF_error = (double **) calloc(m_nenergies, sizeof(double *));
  for (int i = 0; i < m_nenergies; ++i) {
    m_PSF[i] = (double *) calloc(m_nangles, sizeof(double));
    m_PSF_error[i] = (double *) calloc(m_nangles, sizeof(double));
  }
  for (int i = 0; i < m_nenergies; ++i) {
    for (int j = 0; j < m_nangles; ++j) {
      m_PSF[i][j] = intensity[i][j] * m_exposure / n_good_events;
      m_PSF_error[i][j] = intensity_error[i][j] * m_exposure / n_good_events;
    }
  }
  if (m_debug) {
    for (int i = 0; i < m_nenergies; ++i) {
      for (int j = 0; j < m_nangles; ++j) {
  	std::cout << "m_PSF[" << i << "][" << j << "] = " << m_PSF[i][j] <<
 	  std::endl;
  	std::cout << "m_PSF_error[" << i << "][" << j << "] = " <<
 	  m_PSF_error[i][j] << std::endl;
      }
    }
  }

  //---------------------------------------------------------------------------

  // Write the PSFs to a FITS file.
  writeFitsFile();

}

void GtobspsfApp::writeFitsFile() {

  // Get the group of parameters for this application.
  st_app::AppParGroup & pars(getParGroup("gtobspsf"));

  // Fetch the path to the output FITS file.
  std::string output_file = pars["outfile"];
  if (m_debug) {
    std::cout << "output_file = " << output_file << std::endl;
  }

  // Create the output file.
  tip::IFileSvc::instance().createFile(output_file);

  //---------------------------------------------------------------------------

  // Add the PSF extension.

  // Fetch the name of the PSF extension.
  std::string outtable = pars["outtable"];
  if (m_debug) {
    std::cout << "outtable = " << outtable << std::endl;
  }

  // Create the PSF extension in the output file.
  tip::IFileSvc::instance().appendTable(output_file, outtable);

  // Open the PSF table for editing.
  std::auto_ptr<tip::Table>
    psf_table(tip::IFileSvc::instance().editTable(output_file, outtable));

  // Create the Energy field.
  psf_table->appendField("Energy", "1D");

  // Create the Exposure field.
  psf_table->appendField("Exposure", "1D");

  // Create the Psf field.
  std::ostringstream format;
  format << m_nangles << "D";
  psf_table->appendField("Psf", format.str());

  // Create the Psf_err field.
  psf_table->appendField("Psf_err", format.str());

  // Set the table length to the number of energy bins.
  psf_table->setNumRecords(m_nenergies);

  // Fill in each row of the PSF table.
  tip::Table::Iterator row = psf_table->begin();
  tip::Table::Record & record = *row;
  for (int i = 0; i < m_nenergies; ++i) {
    record["Energy"].set(m_bin_energies[i]);
    record["Exposure"].set(m_exposure);
    tip::Table::Vector<double> Psf = record["Psf"];
    tip::Table::Vector<double> Psf_err = record["Psf_err"];
    for (int j = 0; j < m_nangles; ++j) {
      Psf[j] = m_PSF[i][j];
      Psf_err[j] = m_PSF_error[i][j];
    }
    ++row;
  }

  //---------------------------------------------------------------------------

  // Add the bin angle extension.

  // Set the name of the bin angle extension.
  std::string extension_name = "THETA";

  // Create the bin angle extension in the output file.
  tip::IFileSvc::instance().appendTable(output_file, extension_name);

  // Open the bin angle table for editing.
  std::auto_ptr<tip::Table> 
    bin_angle_table(tip::IFileSvc::instance().editTable(output_file,
							extension_name));

  // Create the Theta field.
  bin_angle_table->appendField("Theta", "1D");

  // Set the table length to the number of angle bins.
  bin_angle_table->setNumRecords(m_nangles);

  // Fill in each row of the bin angle table.
  row = bin_angle_table->begin();
  record = *row;
  for (int j = 0; j < m_nangles; ++j) {
    record["Theta"].set(m_bin_angles[j]);
    ++row;
  }

}

//*****************************************************************************

// Create a factory which st_app will use to make the application
// object.
st_app::StAppFactory<GtobspsfApp> g_factory("gtobspsf");
