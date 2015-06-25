// Flux test program that generates a ROOT macro to plot the flux
//

// $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/flux/rootplot.h,v 1.4 2005/06/23 22:52:15 burnett Exp $

// Original author: Theodore Hierath

/**
Test program for graphing the spectrums available through the flux
package.
*/

#include "flux/FluxMgr.h"
#include "flux/FluxSource.h"
#include "flux/SpectrumFactoryTable.h"

#include <fstream>
#include <string>

class rootplot
{
public:
    /// ctor
    rootplot(std::vector<std::string> argv, FluxMgr*);

    rootplot(int argc, char* argv[]);

    const int NUM_BINS;
    const int LOOP;
    
    const double TIME;
    
    const double ENERGY_MIN;
    const double ENERGY_MAX;
    
    void help() {
        std::cout << 
            "   Simple test program for a particle source.\n"
            "   Command line args are \n"
            "      '-sum' sums up the histograms'\n"
            "      '-bins <number of bins for histogram>'\n"
            "      '-events <number of events to create>'\n"
            "      '-min' or '-energy_min' sets minimum energy\n"
            "      '-max' or '-energy_max' sets maximum energy\n"
            "      '-list' lists the available spectra\n"
            "      '-file <file name>' writes energy vs. flux to filename.txt\n"
            "      '-trueflux' graphs flux per steradian\n"
            "      '-no_integrate' is the same as '-trueflux'\n"
            "      '-flux' graphs the flux vs. E instead of E*flux vs E\n"
            "      '-flux_min' manually set lower limit for graph\n"
            "      '-flux_max' manually set upper limit for graph\n"
            "      '-graph <log | semilogx | semilogy | linear>'\n"
            "      '-longsrc <sourcename>' for long-term energy averaging\n"
            "      '-time <time in seconds>' for the flux at time\n"
            "      '-stationary' keeps the satellite from moving (doesn't work for sources\n"
            "                    that use interval to determine the flux)\n"
            "      '-help' for this help"
            << std::endl;
    }
    
    void listSources(const std::list<std::string>& source_list ) {
        std::cout << "List of available sources:" << std::endl;
        for( std::list<std::string>::const_iterator it = source_list.begin(); 
        it != source_list.end(); ++it) { 
            std::cout << '\t'<< *it << std::endl;
        }
    }
    
    void listSpectra() {
        std::cout << "List of loaded Spectrum objects: " << std::endl;
        std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());
        for( std::list<std::string>::const_iterator it = spectra.begin(); 
        it != spectra.end(); ++it) { 
            std::cout << '\t'<< *it << std::endl;
        }
    }
    
      
private:
    void init(std::vector<std::string> argv);

    FluxMgr* m_fm;

};
