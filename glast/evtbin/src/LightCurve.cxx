/** \file LightCurve.cxx
    \brief Encapsulation of a Light curve, with methods to read/write using tip.
    \author James Peachey, HEASARC
*/
#include <memory>
#include <string>

#include "evtbin/Binner.h"
#include "evtbin/LightCurve.h"

#include "st_facilities/Env.h"

#include "facilities/commonUtilities.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

namespace evtbin {

  LightCurve::LightCurve(const std::string & event_file, const std::string & event_table, const std::string & sc_file,
    const std::string & sc_table, const Binner & binner, const Gti & gti): DataProduct(event_file, event_table, gti),
    m_hist(binner) {
    m_hist_ptr = &m_hist;

    // Collect any/all needed keywords from the primary extension.
    harvestKeywords(m_event_file_cont);

    // Collect any/all needed keywords from the GTI extension.
    // But do not fail if GTI isn't there. This is for GBM headers.
    try {
      harvestKeywords(m_event_file_cont, "GTI");
    } catch (...){}

    // Collect any/all needed keywords from the events extension.
    harvestKeywords(m_event_file_cont, m_event_table);

    // Adjust the GTI based on binning information.
    adjustGti(&binner);

    // Update tstart/tstop etc.
    adjustTimeKeywords(sc_file, sc_table, &binner);
  }

  LightCurve::~LightCurve() throw() {}

  void LightCurve::writeOutput(const std::string & creator, const std::string & out_file) const {
    // Standard file creation from base class.
    createFile(creator, out_file, facilities::commonUtilities::joinPath(m_data_dir, "LatLightCurveTemplate"));

    // Open RATE extension of output light curve file. Use an auto_ptr so that the table object
    // will for sure be deleted, even if an exception is thrown.
    std::auto_ptr<tip::Table> output_table(tip::IFileSvc::instance().editTable(out_file, "RATE"));

    // Write DSS keywords to preserve cut information.
    writeDssKeywords(output_table->getHeader());

    // Write the history that came from the events extension.
    writeHistory(*output_table, "EVENTS");

    // The binner from the histogram will be used below.
    const Binner * binner = m_hist.getBinners().at(0);

    // Resize table: number of records in light curve must == the number of bins in the binner.
    output_table->setNumRecords(binner->getNumBins());

    // Need output table iterator.
    tip::Table::Iterator table_itor = output_table->begin();

    // Iterate over bin number and output table iterator, writing fields in order.
    double total_counts=0;
    double total_error_channel=0;
    for (long index = 0; index != binner->getNumBins(); ++index, ++table_itor) {
      // Midpoint time of each bin, from the binner.
      (*table_itor)["TIME"].set(binner->getInterval(index).midpoint());

      // Width of each bin, from the binner.
      (*table_itor)["TIMEDEL"].set(binner->getBinWidth(index));

      // Number of counts in each bin, from the histogram.
      (*table_itor)["COUNTS"].set(m_hist[index]);

      //Keep a running total of binned counts.
      total_counts+=m_hist[index];

      //Statistical Error
      (*table_itor)["ERROR"].set(calcStatErr(m_hist[index]));
    }

    //Check for and if needed make gbm specific correction for deadtime.
    gbmExposure(total_counts, total_error_channel, out_file);

    // Write GTI extension.
    writeGti(out_file);
  }

}
