/**
 * @file IrfHdus.cxx
 * @brief Container for the HDU numbers associated with specified
 * IRF components.
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/irfs/irfUtil/src/IrfHdus.cxx,v 1.1.2.1 2015/02/04 02:01:40 jasercio Exp $
 */

#include <sstream>
#include <stdexcept>

#include "irfUtil/HdCaldb.h"
#include "irfUtil/Util.h"
#include "irfUtil/IrfHdus.h"

namespace irfUtil {

const char * aeff_cnames[3] = {"EFF_AREA", "PHI_DEP", "EFFICIENCY_PARS"};
const char * psf_cnames[3] = {"RPSF", "PSF_SCALING", "FISHEYE_CORR"};
const char * edisp_cnames[2] = {"EDISP", "EDISP_SCALING"};

std::vector<std::string> IrfHdus::s_aeff_cnames(aeff_cnames, aeff_cnames+3);
                                                
std::vector<std::string> IrfHdus::s_psf_cnames(psf_cnames, psf_cnames+2);
                                               
std::vector<std::string> IrfHdus::s_edisp_cnames(edisp_cnames, edisp_cnames+2);

IrfHdus::IrfHdus(const std::string & irf_name,
                 const std::string & event_type,
                 const std::vector<std::string> & cnames) 
   : m_cnames(cnames), m_convType(0) {
   if (event_type == "BACK") {
      /// The m_convType value is required by the pre-Pass 8 PSF
      /// classes to determine which PSF scaling parameters to use.
      /// It need only differ from zero for BACK section IRFs, but for
      /// backwards-compatibility, Pass 8 and later IRFs must support
      /// m_convType=1.
      m_convType = 1;
   }
   irfUtil::HdCaldb hdcaldb("GLAST", "LAT");

   for (size_t i(0); i < cnames.size(); i++) {
      std::vector<std::string> filenames;
      std::vector<int> hdus;
      hdcaldb.getFiles(filenames, hdus, event_type, cnames[i], irf_name);
      FilenameHduPairs_t fh_pairs;
      for (size_t j(0); j < hdus.size(); j++) {
         std::ostringstream extname;
         extname << hdus[j];
         fh_pairs.push_back(std::make_pair(filenames.at(j), extname.str()));
      }
      m_file_hdus[cnames[i]] = fh_pairs;
   }

   std::map<std::string, std::pair<unsigned int, std::string> > evtype_mapping;
   irfUtil::Util::get_event_type_mapping(irf_name, evtype_mapping);
   m_bitPos = evtype_mapping[event_type].first;
}

typedef std::vector< std::pair<std::string, std::string> > FilenameHduPairs_t;

const FilenameHduPairs_t & 
IrfHdus::operator()(const std::string & cname) const {
   const std::map<std::string, FilenameHduPairs_t>::const_iterator 
      it(m_file_hdus.find(cname));
   if (it == m_file_hdus.end()) {
      throw std::runtime_error("latResponse::IrfHdus: " + cname 
                               + " not found.");
   }
   return it->second;
}

unsigned int IrfHdus::bitPos() const {
   return m_bitPos;
}

int IrfHdus::convType() const {
   return m_convType;
}

size_t IrfHdus::numEpochs() const {
   return m_file_hdus.begin()->second.size();
}

const std::vector<std::string> & IrfHdus::cnames() const {
   return m_cnames;
}

IrfHdus IrfHdus::aeff(const std::string & irf_name,
                      const std::string & event_type) {
   return IrfHdus(irf_name, event_type, s_aeff_cnames);
}

IrfHdus IrfHdus::psf(const std::string & irf_name,
                     const std::string & event_type) {
   return IrfHdus(irf_name, event_type, s_psf_cnames);
}

IrfHdus IrfHdus::edisp(const std::string & irf_name,
                       const std::string & event_type) {
   return IrfHdus(irf_name, event_type, s_edisp_cnames);
}

} // namespace latResponse
