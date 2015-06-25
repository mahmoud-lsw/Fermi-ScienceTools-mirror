/** @file rootplot.h

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/src/rootplot/rootplot.cxx,v 1.15 2011/05/20 16:14:56 heather Exp $
*/
#include "flux/rootplot.h"

#include "rootEnergyHist.h"
#include "rootAngleHist.h"

#include "CLHEP/Geometry/Vector3D.h"
using astro::GPS;

namespace {

    class Maglat {
    public:
        Maglat():m_mean(0), m_min(99), m_max(-99), m_count(0){}
        void accum(double x){
            if(x< m_min) m_min=x;
            if(x>m_max) m_max =x ;
            m_mean += x;
            m_count++;
        }
        void dump(){
            std::cout << "\nMagnetic latitude:"
                << "\n\tmin\t" << m_min
                << "\n\tmax\t" << m_max
                << "\n\tmean\t"<< m_mean/m_count
                << std::endl;
        }
    private:
        double m_mean, m_min, m_max;
        int m_count;
    };

}

rootplot::rootplot(int argc, char* argv[])
: NUM_BINS(30),LOOP(30000),
TIME(0.01), 
ENERGY_MIN(0.01*1000.), 
ENERGY_MAX(100.0*1000.), m_fm(0)
{
    std::vector<std::string> args;
    for( int i =1; i< argc; ++i) 
        args.push_back(argv[i]);
    std::vector<std::string > test;
    FluxMgr fm(test);
    m_fm = &fm;
    init(args);

}


rootplot::rootplot(std::vector<std::string> argv, FluxMgr* fm)
: NUM_BINS(30),LOOP(30000),
TIME(0.0), ENERGY_MIN(0.01*1000.), ENERGY_MAX(100.0*1000.)
,m_fm(fm)
{
    init(argv);
}

void rootplot::init(std::vector<std::string> argv)
{

    int argc = argv.size();
    static std::string default_arg="default";
    static std::string default_graph="log";

    int num_bins = NUM_BINS;  // Initialize to default number of bins
    int loop = LOOP; // Initialize to default number of iterations
    int num_sources = 0;
    int num_longsources = 0;
    int current_arg = 0;
    int longtime=1;
    int longtimemax=20;
    double time=TIME;  //time to use for flux and rate functions
    double energy_min = ENERGY_MIN;
    double energy_max = ENERGY_MAX;
    double expand(1e4); // expansion factor?
    bool use_integrated_flux = true;
    bool use_flux = false;
    bool use_flux_min = false;
    bool use_flux_max = false;
    bool write_to_file = false;
    bool sum_mode = false;
    bool longterm=false;
    bool stationary=false;
    double flux_min;
    double flux_max;
    std::string arg_name(default_arg);
    std::string output_file_name;
    std::string ylabel_flux   = "Flux (particles/m^2/s/MeV/sr)";
    std::string ylabel_eflux  = "E*Flux (particles/m^2/s/sr)";
    std::string ylabel_integ_flux  = "Flux (particles/m^2/s/MeV)"; // Flux integrated over solid angle
    std::string ylabel_integ_eflux = "E*Flux (particles/m^2/s)";   // E*Flux integrated over solid angle

    std::vector<std::string> sources;
    std::vector<int> longsources;

    //    flux_load();

    //FluxMgr fm(sources); 
    FluxMgr & fm = *m_fm;


    // Process Command Line Arguments

    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "flux test program: type 'rootplot -help' for help" << std::endl;
    std::cout << ( ( argc == 0)?  " No command line args, using defaults"
        :  "") << std::endl;



    while(current_arg < argc)
    {
        arg_name = argv[current_arg];
        if("-help" == arg_name || "help" == arg_name) 
        { 
            help(); 
            return;// 0; 
        }
        else if("-bins" == arg_name) 
            num_bins = atoi(argv[++current_arg].c_str());
        else if("-events" == arg_name) 
            loop = atoi(argv[++current_arg].c_str());
        else if("-min" == arg_name || "-energy_min" == arg_name) 
            energy_min = atof(argv[++current_arg].c_str());
        else if("-max" == arg_name || "-energy_max" == arg_name) 
            energy_max = atof(argv[++current_arg].c_str());
        else if("-flux_min" == arg_name)
        {
            use_flux_min = true;
            flux_min = atof(argv[++current_arg].c_str());
        }
        else if("-flux_max" == arg_name)
        {
            use_flux_max = true;
            flux_max = atof(argv[++current_arg].c_str());
        }
        else if("-flux" == arg_name)
            use_flux = true;
        else if("-trueflux" == arg_name || "-no_integrate" == arg_name)
            use_integrated_flux = false;
        else if("-file" == arg_name)
        {
            write_to_file = true;
            output_file_name = argv[++current_arg];
        }
        else if("-sum" == arg_name)
            sum_mode = true;
        else if("-list" == arg_name) 
        { 
            listSources(fm.sourceList());
            listSpectra(); 
            return; 
        }
        else if("-graph" == arg_name) 
            default_graph = argv[++current_arg];
        else if("-liny" == arg_name) default_graph = "semilogx";
        else if("-longsrc" == arg_name) {
            sources.push_back(argv[++current_arg]);
            //put the number of this source into the list for reference
            longsources.push_back(num_sources);
            num_sources++;
        }
        else if("-time" == arg_name) {
            time = atof(argv[++current_arg].c_str());
            std::cout<<" TIME = "<< time << std::endl;
        }
        else if("-stationary" == arg_name) {
            stationary = true;
        }
        else if('-' == arg_name[0]) {std::cerr << "Unrecognized option "<< arg_name << ", -help for help" << std::endl;}
        else
        {
            sources.push_back(arg_name);
            num_sources++;
        }
        current_arg++;
    }


    // Use default source if no source was specified in arguments
    if(0 == sources.size())
    {
        sources.push_back(default_arg);
        num_sources++;
    }

    rootEnergyHist energy_hist(num_bins,energy_min,energy_max, "rootplot.cxx");
    rootAngleHist angle_hist(num_bins);

    // Process all the sources
    for(int i = 0; i < num_sources; i++)
    {
        int j;

        //decide whether or not this run should be longterm
        longterm = false;
        std::vector<int>::iterator longiter;
        for( longiter=longsources.begin(); longiter!=longsources.end() ;longiter++){
            if(*longiter==i) longterm=true;
        }

        if(longterm)
            fm.setExpansion(-1.);


        if((false == sum_mode && false==longterm)||(true==longterm && (longtime==1))) 
        {
            // Reset for new source
            energy_hist.reset();
            angle_hist.reset();
        }

        // Make sure positions are initialized
        GPS::instance()->synch();
        fm.pass(0.);

        std::string sourcename(sources[i])
                  , sourcedisplayname(sources[i]);
        size_t colon (sourcename.find(":"));
        if(colon!= std::string::npos){
            std::string check1 = sourcename.substr(0, colon);
            std::string check2 = sourcename.substr(colon+1);
            sourcedisplayname= check2;
            sourcename = check1;
        }
        EventSource *e = fm.source(sourcename);
#if 0

        if(longterm){
            fm.pass(2.);
            time+=2.;
        }else{time=0;}

#endif

        if( 0==e ) {std::cerr << "Source \"" << sources[i] << "\" not found: -list for a list" << std::endl;
        return;}

        energy_hist.setGraphType(default_graph.c_str());
        energy_hist.setTitle( sourcedisplayname );

        energy_hist.setXLabel("Kinetic Energy (MeV)");

        if(true == use_flux)
        {
            energy_hist.setFluxMode();
            if(true == use_integrated_flux)
                energy_hist.setYLabel(ylabel_integ_flux);
            else
                energy_hist.setYLabel(ylabel_flux);
        }
        else
        {
            if(true == use_integrated_flux)
                energy_hist.setYLabel(ylabel_integ_eflux);
            else
                energy_hist.setYLabel(ylabel_eflux);
        }

        if(true == use_flux_min)
            energy_hist.setFluxMin(flux_min);

        if(true == use_flux_max)
            energy_hist.setFluxMax(flux_max);

        angle_hist.setGraphType(default_graph.c_str());
        angle_hist.setTitle( sourcedisplayname );
        angle_hist.setPhiXLabel("Angle (degrees)");
        angle_hist.setPhiYLabel("Particles");
        angle_hist.setThetaXLabel("Cos(Theta)");
        angle_hist.setThetaYLabel("Particles");


        std::cout << sources[i] << std::endl;

        GPS* gps = GPS::instance();
        gps->time(time);
        std::pair<double,double> loc=fm.location();
        std::cout << "Lat/Lon:  " << loc.first << "   " << loc.second << std::endl;

        //	  std::cout << "orbit angle=" << GPS::instance()-> << "orbit phase=" << << std::endl;

        std::cout << "Initial (p/s/m^2): " << e->rate(time)/e->totalArea() << std::endl;
        std::cout << "Initial (p/s/m^2/sr): " << e->flux(time) << std::endl;

        double time2 = 0;  // time used for scaling the graphs
        Maglat magstat;
        for(j = 0; j < loop; j++) 
        {
            EventSource *f = e->event(time2); //time);
            //increment the time
            double timeadd = e->interval();
            time2 += timeadd;

            if(!stationary)
            {
                time += timeadd* expand;
                fm.pass(timeadd* expand);
                //fm.time(time);
            }
            gps->time(time); // set to expandec time
            gps->notifyObservers();
            double geolat = gps->earthpos().geolat();
            magstat.accum(geolat);
 
//            double cutOffRigidity = gps->expansion();
            HepGeom::Vector3D<double> dir = f->launchDir();
            double energy = f->energy();
            gps->time(time2); // reset

            double cos_theta = dir.z();

            double phi = atan2(dir.y(),dir.x());
            if(phi < 0)
                phi = 2*M_PI + phi;

            energy_hist.store(energy);
            angle_hist.storeTheta(cos_theta);
            angle_hist.storePhi(phi);

            if(j % 1000 == 0) {std::cout << "\r" << j << ": " << energy << "...";
            }
        }
        magstat.dump();

        std::cerr << "\n";

        double scale_factor;


        double flux=(loop/time2)/e->totalArea();
        // There's probably a better way to get at the solid angle but
        // e->solidAngle doesn't seem to be initialized.
        double solidangle = e->rate(time)/(e->flux(time)*e->totalArea());

        // These scale factors are approximately correct for sources that have
        // spectrums that change over time.  For static sources these will work 
        // very well.  To do: Apply scale factor to each event individually rather
        // than an average scale factor to them all at the end.
        if(true == use_integrated_flux)
            scale_factor = flux/loop*num_bins/log(energy_max/energy_min);
        else
            scale_factor = flux/loop*num_bins/log(energy_max/energy_min)/solidangle;

        std::cout << "Average (p/s/m^2): " << flux  << std::endl;
        std::cout << "Average (p/s/m^2/sr): " << flux/solidangle << std::endl;


        // Scale factor must be applied before the draw member function can be called
        energy_hist.apply(scale_factor);
        angle_hist.apply(scale_factor);

        for(j = 0; j < num_bins; j++) 
        {
            std::cout << energy_min*pow(10.0,((j + 0.5)/num_bins) * energy_hist.retrieveRange() ) 
                << "   \t" << energy_hist.retrieveFlux(j) << "\t";
            std::cout << std::endl;         
        }


        if(false == sum_mode && false==longterm)
        {
            angle_hist.draw(1,"begin",i,num_sources);
            energy_hist.draw(1,"end",i,num_sources);

            //delete e;
            std::cout << "Normal method" << std::endl;        
        }
        else if(true == sum_mode && i == num_sources - 1)
        {
            angle_hist.draw(1,"begin",0,1);
            energy_hist.draw(1,"end",0,1);

            //delete e;
            std::cout << "Sum Mode method" << std::endl;
        }
#if 0
        else if(longterm==true && longtime<=longtimemax)
        {
            longtime++;
            i--;
            std::cout << "Longterm run number " << longtime << std::endl; 
        }
#endif
        else
        {       
            if(true == write_to_file)
            {
                std::ofstream output_file;
                output_file_name += ".txt";
                output_file.open(output_file_name.c_str());

                energy_hist.writeFile(1./(longtime),output_file);
                output_file.close();
            }
            angle_hist.draw(1./double(longtime),"begin",i,num_sources);
            energy_hist.draw(1./double(longtime),"end",i,num_sources);

            delete e;
            longtime=1; 

        }

        std::cout << std::endl;

    } // for all sources   

   return;
}
