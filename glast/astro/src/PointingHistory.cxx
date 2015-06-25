/** @file PointingHistory.cxx
    @brief implement PointingHistory

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/astro/src/PointingHistory.cxx,v 1.25 2013/06/28 20:44:17 jperkins Exp $

    */

#include "astro/PointingHistory.h"
using namespace astro;

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/Extension.h"

#include <fstream>
#include <sstream>
#include <memory>
#include <cmath>
using CLHEP::Hep3Vector;

namespace {
    static double time_tol(10); // seconds allow beyond the end
    static double max_interval(65); // more than this triggers warning message
//    int max_warning_message_count(100); // suppress messages if more than this  
    int max_warning_message_count(0); // suppress messages if more than this  
}

PointingHistory::PointingHistory(const std::string& filename, double offset)
: m_selected(-1)
, m_startTime(-1)
{
    if( haveFitsFile(filename) ){
        readFitsData(filename);
    }else{
        readTextData(filename, offset);
    }
}
void PointingHistory::readTextData(std::string filename, double offset)
{
    std::ifstream input_file;
    input_file.open(filename.c_str());

    if(!input_file.is_open()) {
        std::cerr << "ERROR:  Unable to open:  " << filename << std::endl;
        throw std::invalid_argument("PointingHistory: could not open pointing history file");
    }else{
        //initialize the key structure:
        while (!input_file.eof()){
            std::string line; std::getline(input_file, line);
            double start; std::stringstream buf(line); buf >>start;
            double x,y,z;
            buf >> x >> y >> z;
            Hep3Vector position(x,y,z);

            double raz, decz, rax, decx;
            buf >> raz >> decz >> rax >> decx;
            Quaternion orientation(SkyDir(raz,decz)(), SkyDir(rax,decx)());

            double razenith, deczenith;
            buf >> razenith >> deczenith; // ignore since redundant with position
            double lat, lon, alt;
            buf >> lon >> lat >> alt;     // these are ignored, too.

            EarthCoordinate earth(position, start);

            m_endTime = start + offset;
            m_data[m_endTime] = PointingInfo(position, orientation, earth);
            if( m_startTime<0) m_startTime=m_endTime;

        }
    }
}

const astro::PointingInfo& PointingHistory::operator()(double time)const throw(TimeRangeError){

    if( time!=m_selected){

        std::map<double,astro::PointingInfo>::const_iterator iter=m_data.upper_bound(time);
        bool tooEarly( time< (*(m_data.begin())).first );
        bool tooLate( iter==m_data.end() && time > m_endTime);

        if(tooEarly && time==0) {
            // special case: just return starting value
            time = (*(m_data.begin())).first;
            tooEarly = false;
        }

	if ( tooEarly || tooLate) {
	  std::ostringstream message;
	  if ( tooEarly ) {
	    message << "PointingHistory: Time out of Range!:\n"
		    << "Time (" << static_cast<int>(time)
		    << ") occurs " << m_startTime - time << "s"
		    << " before the range of times in the pointing database: ("
		    << static_cast<int>(m_startTime) << ", "<< static_cast<int>(m_endTime) << ")"
		    << std::endl;
	  } else if ( tooLate ) {
	    message << "PointingHistory: Time out of Range!:\n"
		    << "Time (" << static_cast<int>(time)
		    << ") occurs " << time - m_endTime << "s"
		    << " after the range of times in the pointing database: ("
		    << static_cast<int>(m_startTime) <<", "<< static_cast<int>(m_endTime) <<")"
		    << std::endl;
	  }
	  throw TimeRangeError(message.str());
        }
        const PointingInfo & h2 = iter->second;
        double time2( iter->first);
        --iter;
        double time1(iter->first);
        if( time-time1 > max_interval && max_warning_message_count>0 ){ 
            --max_warning_message_count;
            std::cerr << "Warning: " 
                << int(time+0.5) << " seems to be in an invalid interval: it is " 
                << (time-time1) << " seconds beyond the start of the current FT2 entry" 
                << std::endl;
        }
        const PointingInfo & h1 =iter->second;
        double prop( (time-time1)/(time2-time1) );

        m_currentPoint = h1.interpolate(h2, prop, time);

    }
    return m_currentPoint;
}

bool PointingHistory::haveFitsFile(std::string filename) const {

     // Determine a type of the input file, FITS or TEXT.
     bool is_fits = true;
     try {
       // Try opening a primary extension of a FITS file.
         std::auto_ptr<const tip::Extension>ext(tip::IFileSvc::instance().readExtension(filename, "0"));
     } catch (const tip::TipException &) {
       is_fits = false;
     }
     return is_fits;
}

void PointingHistory::readFitsData(std::string filename) {

    std::string sctable("SC_DATA"); //?
    const tip::Table * scData = 
        tip::IFileSvc::instance().readTable(filename, sctable);
    double start_time, stop_time, raz, decz, rax, decx, lat,lon;
    std::vector<double> sc_pos(3);
    tip::Table::ConstIterator it = scData->begin();
    tip::ConstTableRecord & interval = *it;
    double maxlatdiff(0);
    double last_start(0), last_stop(0);
    bool have_quaternion(false);
    double qsj_1, qsj_2, qsj_3, qsj_4;
#if 0 // back this out
    try{
        interval["QSJ_1"].get(qsj_1);
        have_quaternion= true;
    }catch(const std::exception&){
    }
#endif
    for ( ; it != scData->end(); ++it) {

        Quaternion orientation;
        if( have_quaternion ) {
            interval["QSJ_1"].get(qsj_1);
            interval["QSJ_2"].get(qsj_2);
            interval["QSJ_3"].get(qsj_3);
            interval["QSJ_4"].get(qsj_4);
            orientation = Quaternion(Hep3Vector(qsj_1,qsj_2,qsj_3), qsj_4);

        }else{
            interval["ra_scz"   ].get(raz);
            interval["dec_scz"  ].get(decz);      
            interval["ra_scx"   ].get(rax);
            interval["dec_scx"  ].get(decx);  
            astro::SkyDir zaxis(raz, decz);
            astro::SkyDir xaxis(rax, decx);
            orientation = Quaternion(zaxis(), xaxis());
        }
        interval["lat_geo"  ].get(lat);
        interval["lon_geo"  ].get(lon);
        interval["start"    ].get(start_time);
        interval["stop"     ].get(stop_time);
        interval["sc_position"].get(sc_pos);
        CLHEP::Hep3Vector position(sc_pos[0]/1e3, sc_pos[1]/1e3, sc_pos[2]/1e3);
        EarthCoordinate earthpos(position, start_time);

        // Read the properties of the LAT for this interval.
        int lat_mode;
        int lat_config;
        int data_qual;
        double livetime;
        double rock_angle;
        bool in_saa;
        interval["lat_mode"].get(lat_mode);
        interval["lat_config"].get(lat_config);
        interval["data_qual"].get(data_qual);
        interval["livetime"].get(livetime);
        interval["rock_angle"].get(rock_angle);
        interval["in_saa"].get(in_saa);

        // check consistency of latitude, longitude: EarthCoordinate computes from the MET and position
        double check_lat(earthpos.latitude()-lat), check_lon(earthpos.longitude()-lon);
        if( check_lat>maxlatdiff) maxlatdiff = check_lat;

        static double lat_tol(0.5), lon_tol(0.5); // was 0.2, 0.01
        static int error_count(50); // allow this many

        if( (fabs(check_lat)>lat_tol || fabs(check_lon)>lon_tol && std::abs(lon)<179) && error_count>0 ){
            std::stringstream error; 
            error << "PointingHistory::readFitsData: apparent inconsistency for Earth position, time=" 
                << start_time
                << ", lat, lon diff: " << check_lat << ", " << check_lon;
            std::cerr << error.str() << std::endl;
#if 0 // disable for now
            throw std::runtime_error(error.str());
#else // allow, but limit number of messages
            error_count--;
            if (error_count==0){
                std::cerr << "\tSuppressing error messages: beware that LAT/LON columns are unreliable" << std::endl;
            }
#endif
        }
        // add an entry for the end of a run -- make identical (a small kluge)
       if( last_stop>0 && start_time != last_stop){
            m_data[last_stop] = m_data[last_start]; 
       }
       m_data[start_time] = 
           PointingInfo(position, orientation, earthpos, 
                        LatProperties(lat_mode, lat_config, data_qual,
                                      livetime, start_time, stop_time,
                                      rock_angle, in_saa));
    
       if( m_startTime<0) m_startTime = start_time;
       last_start = start_time;
       last_stop = stop_time;
    }
    m_endTime = stop_time+time_tol;
    // an extra entry to allow query for the last interval +slop
    m_data[m_endTime] = m_data[start_time];
}
