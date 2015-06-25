/**
 * @file SpectralTransient.h
 * @brief A flaring source whose light curve and spectral variability
 * is given by a template file.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/genericSources/genericSources/SpectralTransient.h,v 1.12 2007/12/27 06:32:48 jchiang Exp $
 */

#ifndef genericSources_SpectralTransient_h
#define genericSources_SpectralTransient_h

#include "flux/Spectrum.h"

namespace IRB {
   class EblAtten;
}

/**
 * @class SpectralTransient
 *
 * @brief A flaring source whose light curve shape is given by a
 * template file.  The duration and mean photon flux are given as
 * parameters.  The spectrum during each interval defined in the
 * template file is given as a broken power-law.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/celestialSources/genericSources/genericSources/SpectralTransient.h,v 1.12 2007/12/27 06:32:48 jchiang Exp $
 */

class SpectralTransient : public Spectrum {

public:

   SpectralTransient(const std::string & params);
   
   virtual ~SpectralTransient();

   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi) const {
      (void)(xi);
      return m_currentEnergy;
   }

   /// @return Particle type, "gamma" is the default.
   virtual const char * particleName() const {
      return m_particle.c_str();
   }

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "SpectralTransient";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time) {
      (void)(time);
      return m_currentEnergy;
   }

private:

   double m_flux;
   double m_tstart;
   double m_tstop;
   double m_emin;
   double m_emax;
   int m_lc;
   float m_z;
   int m_logParabola;
   bool m_specFile;

   IRB::EblAtten * m_tau;

   double m_tauScale;

   std::string m_particle;

   double m_currentEnergy;

   void readModel(std::string templateFile);

   void readSpectrum(std::string specFile);

   std::vector<double> m_energies;
   std::vector<double> m_dnde;
   std::vector<double> m_cumulativeDist;

   double drawEnergy() const;

   class ModelInterval {
   public:

      ModelInterval() {}

      /// Read data members from a line in the template file.
      ModelInterval(const std::string & line, double emin, double emax,
                    int useLogParabola=0);

      /// Pass the data members via an ordered vector.
      ModelInterval(const std::vector<double> & data,
                    double emin, double emax, int useLogParabola=0);

      /// Fractional start time of the interval; the entire light curve 
      /// will be rescaled to fit the interval [m_tstart, m_tstop]
      double startTime;

      /// Fractional stop time.
      double stopTime;

      /// Energy integrated photon flux (#/cm^2-s)
      double flux;

      /// Lower energy photon index
      double gamma1;

      /// Higher energy photon index
      double gamma2;

      /// Break energy (MeV)
      double ebreak;

      /// Draw a photon energy (MeV)
      double drawEnergy(double emin, double emax) const;

      void fillCumulativeDist(double emin, double emax);

      void clearCumulativeDist() {
         m_integral.clear();
      }
   private:
      double m_lowerFraction;
      void brokenPowerLawFractions(double emin, double emax);
      bool drawBelowBreak(double xi) const {
         return xi < m_lowerFraction;
      }

      static std::vector<double> s_energies;
      std::vector<double> m_integral;

      double logParabola(double energy) const;
      void fillEnergies(double emin, double emax);
   };

   std::vector<ModelInterval> m_lightCurve;

   std::vector< std::pair<double, double> > m_eventCache;

//   std::vector<ModelInterval>::const_iterator m_currentInterval;
   std::vector<ModelInterval>::iterator m_currentInterval;

   void readLightCurve(const std::string & templateFile);
   void readFitsLightCurve(const std::string & templateFile);

   void rescaleLightCurve();

   void fillEventCache(double time);

   double drawEnergy(const ModelInterval & interval) const;

   static bool compareEventTime(const std::pair<double, double> & x,
                                const std::pair<double, double> & y) {
      return x.first > y.first;
   }
};

#endif // genericSources_SpectralTransient_h
