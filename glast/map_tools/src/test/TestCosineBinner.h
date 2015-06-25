
/** @file TestCosineBinner.h
@brief test class for CosineBinner

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/map_tools/src/test/TestCosineBinner.h,v 1.3 2007/12/11 05:06:58 burnett Exp $


*/
#include "healpix/CosineBinner.h"
#include <stdexcept>


class TestCosineBinner {
public:
    // simple function to test the integral

    class CosFun{
    public:
        CosFun(double cmin=healpix::CosineBinner::cosmin()):m_cmin(cmin){}
        double operator()(double cth)const{
            return (cth-m_cmin)/(1-m_cmin);
        }
        double m_cmin;
    };

    TestCosineBinner(std::ostream& out= std::cout)
    {
        using healpix::CosineBinner;

        // test the cos binner
        out << "\nTesting CosineBinner: " << std::endl;

        CosineBinner::setBinning(-1., 10, false);
        out << "Linear bins from 1.0 to "<< CosineBinner::cosmin()<< std::endl;
        CosineBinner linear;
        for( double z=CosineBinner::cosmin()+0.005; z< 1.0; z+=0.01) linear[z] +=0.1f; // make uniform

        for( CosineBinner::const_iterator bit=linear.begin(); bit!= linear.end();++bit){
            out << linear.costheta(bit) <<  " " << *bit << std::endl;
        }
        // and a function to integrate over
        double alinear = linear(CosFun());

        out << "sqrt-weighting bins from 1.0 to "<< CosineBinner::cosmin()<< std::endl;
        CosineBinner::setBinning(-1., 10, true);
        CosineBinner sqrtwt;
        for( double z=CosineBinner::cosmin()+0.005; z< 1.0; z+=0.01) sqrtwt[z] +=0.1f; // make uniform

        for(CosineBinner::const_iterator bit=sqrtwt.begin() ; bit!= sqrtwt.end();++bit){
            out << sqrtwt.costheta(bit) <<  " " << *bit << std::endl;
        }

        double asqrt = sqrtwt(CosFun());
        out << "integral with CosFun, linear: " << alinear <<", sqrt: " << asqrt << std::endl;
        if (fabs(alinear - asqrt)> 0.1) 
            throw std::runtime_error("TestCosineBinner integrals not close enough: see output");
    }
};
