/** @file Healpix.cxx
    @brief Healpix class implementation with code from WMAP

    @author B. Lesnick 
    $Header: /glast/ScienceTools/glast/healpix/src/Healpix.cxx,v 1.1.1.4.6.5 2015/04/26 16:11:50 jasercio Exp $
*/
/* Local Includes */

#include "healpix/Healpix.h"
#include "healpix/base/healpix_base.h"
#include "healpix/base/pointing.h" ///< NASA library for ra,d
#include "healpix/base/arr.h"

///< NASA healpix class
/* Standard Includes */
#include <numeric> // for accumulate
#include <stdexcept>

// anonymous namespace for helper functions
namespace {
    inline double radians(double x){return x*M_PI/180.;}
    inline double degrees(double x){return x*180./M_PI;}
}

using namespace healpix;


//=========================================================================================
//  C++ interface funtions
//=========================================================================================

Healpix::Healpix(long nside, healpix::Healpix::Ordering ord, astro::SkyDir::CoordSystem coordsys)
    : m_coordsys(coordsys)
    , m_nside(nside)
    , m_ord(ord)
{
}
void Healpix::pix2ang(long index, double &theta, double &phi)const
{
    Healpix_Base hp(m_nside,static_cast<Healpix_Ordering_Scheme>(m_ord),SET_NSIDE); 
    pointing point = hp.pix2ang(index);
    theta = point.theta;
    phi = point.phi;
}

void Healpix::ang2pix(double theta, double phi, long &index)const
{
    Healpix_Base hp(m_nside,static_cast<Healpix_Ordering_Scheme>(m_ord),SET_NSIDE); 
    index = hp.ang2pix(pointing(theta,phi));
}

Healpix::Pixel::Pixel(const astro::SkyDir &dir, const Healpix& hp)
: m_healpix(&hp)
{
    // get theta, phi (radians) in appropriate coordinate system
	double theta(M_PI/2), phi(M_PI/180);
    if( hp.coordsys()==astro::SkyDir::EQUATORIAL){
        theta -= radians(dir.dec());
        phi   *= dir.ra();
    }else{  // Galactic
        theta -= radians(dir.b());
        phi   *= dir.l();
    }
    // and look up the pixel number
    m_healpix->ang2pix(theta, phi,m_index);
}

Healpix::Pixel::operator astro::SkyDir ()const
{
    double theta,phi;
    m_healpix->pix2ang( m_index,theta,phi);
    // convert to ra, dec (or l,b)
    return astro::SkyDir( degrees(phi), degrees(M_PI/2-theta), m_healpix->coordsys() );
}

void Healpix::Pixel::neighbors(std::vector<Healpix::Pixel> & p) const
{
    fix_arr<int,8> result;
    
    p.clear();
#if 0 // now only slower?
    if (!(this->m_healpix->nested()))
        throw std::runtime_error("Nested ordering required to determine neighbors.");
#endif
    const Healpix &hpx = *this->m_healpix;
    Healpix_Base hpb(hpx.nside(),static_cast<Healpix_Ordering_Scheme>(hpx.ord()),SET_NSIDE); 
    hpb.neighbors(m_index, result);
    if(result[result.size()-1]  != -1) {
            p.push_back(Healpix::Pixel(result[result.size()-1], *(this->m_healpix)));
    }
    for (int i = 0; i < result.size()-1; ++i)
    {
        if(result[i]  != -1) {
            p.push_back(Healpix::Pixel(result[i], *(this->m_healpix)));
        }
    }   
}

double Healpix::integrate(const astro::SkyFunction& f)const
{
    return std::accumulate(begin(), end(), 0., Integrand(f));
}

void Healpix::findNeighbors(int index, std::vector<int> &p) const
{
   
#if 0  // now only slower?
    if (!(nested()))
        throw std::invalid_argument("Healpix::findNeighbors -- Nested ordering required to determine neighbors.");
#endif
    p.clear();

    fix_arr<int,8> result;// local copy of the list
    Healpix_Base hpb(nside(),static_cast<Healpix_Ordering_Scheme>(ord()),SET_NSIDE); 
    hpb.neighbors(index, result);

    long n[8];
    int nit=0;
    if(result[result.size()-1]!=-1) {
            n[0]=result[result.size()-1];
            nit++;
    }
    for (int i=0;i<result.size()-1;i++) {
        if(result[i]!=-1) {
            n[nit]=result[i];
            nit++;
        }
    }
    // now insert this list into the output vector
    std::copy(n,n+nit, std::back_insert_iterator<std::vector<int> >(p));
}

 ///@brief the number of sides 
long Healpix::nside()const{return m_nside; }
    ///@brief the number of pixels
long Healpix::npix()const{return 12*nside()*nside();}

    ///@brief the number of pixels, as the size.
size_t Healpix::size()const{return npix();}

Healpix::Ordering Healpix::ord()const{return m_ord;}

bool Healpix::nested()const{return ord()==NESTED;}

void Healpix::query_disc (const astro::SkyDir dir, double radius, std::vector<int> & v) const
{
    v.clear();
	double theta(M_PI/2), phi(M_PI/180);
    if( m_coordsys==astro::SkyDir::EQUATORIAL){
        theta -= radians(dir.dec());
        phi   *= dir.ra();
    }else{  // Galactic
        theta -= radians(dir.b());
        phi   *= dir.l();
    }
    Healpix_Base hpb(this->nside(),static_cast<Healpix_Ordering_Scheme>(this->ord()),SET_NSIDE); 
	hpb.query_disc(pointing(theta, phi), radius, v);
}
