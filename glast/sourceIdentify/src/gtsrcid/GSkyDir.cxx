/*------------------------------------------------------------------------------
Id ........: $Id: GSkyDir.cxx,v 1.1.1.1 2010/10/01 17:36:15 elwinter Exp $
Author ....: $Author: elwinter $
Revision ..: $Revision: 1.1.1.1 $
Date ......: $Date: 2010/10/01 17:36:15 $
--------------------------------------------------------------------------------
$Log: GSkyDir.cxx,v $
Revision 1.1.1.1  2010/10/01 17:36:15  elwinter
Import of ScienceTools-v9r18p4-slac-20101001 from SLAC

Revision 1.1  2010/04/16 16:16:19  jurgen
Implement HEALPix interface to read counterpart density maps

------------------------------------------------------------------------------*/
/**
 * @file GSkyDir.cxx
 * @brief Implements GSkyDir methods.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <cmath>
#include "GSkyDir.h"

/* Namespace definition __________________________________________________ */
namespace sourceIdentify {

/* __ Prototype __________________________________________________________ */
double modulo(double v1, double v2);
void   euler(const int& type, const double& xin, const double &yin, 
             double* xout, double *yout);

/* __ Constants __________________________________________________________ */
const double pi         = 3.1415926535897931159979635;
const double twopi      = 2.0 * pi;
const double fourpi     = 4.0 * pi;
const double deg2rad    =  0.0174532925199432954743717;
const double rad2deg    = 57.295779513082322864647722;


/*==========================================================================
 =                                                                         =
 =                      GSkyDir constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GSkyDir::GSkyDir()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Sky direction from which class should be instantiated.
 ***************************************************************************/
GSkyDir::GSkyDir(const GSkyDir& dir)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyDir::~GSkyDir()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GSkyDir operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 ***************************************************************************/
GSkyDir& GSkyDir::operator= (const GSkyDir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(dir);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GSkyDir public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set equatorial sky direction (radians)
 *
 * @param[in] ra Right Ascension in radians.
 * @param[in] dec Declination in radians.
 ***************************************************************************/
void GSkyDir::radec(const double& ra, const double& dec)
{
    // Set attributes
    m_has_lb    = 0;
    m_has_radec = 1;

    // Set direction
    m_ra  = ra;
    m_dec = dec;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set equatorial sky direction (degrees)
 *
 * @param[in] ra Right Ascension in degrees.
 * @param[in] dec Declination in degrees.
 ***************************************************************************/
void GSkyDir::radec_deg(const double& ra, const double& dec)
{
    // Set attributes
    m_has_lb    = 0;
    m_has_radec = 1;

    // Set direction
    m_ra  = ra  * deg2rad;
    m_dec = dec * deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set galactic sky direction (radians)
 *
 * @param[in] l Galactic longitude in radians.
 * @param[in] b Galactic latitude in radians.
 ***************************************************************************/
void GSkyDir::lb(const double& l, const double& b)
{
    // Set attributes
    m_has_lb    = 1;
    m_has_radec = 0;

    // Set direction
    m_l = l;
    m_b = b;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set galactic sky direction (degrees)
 *
 * @param[in] l Galactic longitude in degrees.
 * @param[in] b Galactic latitude in degrees.
 ***************************************************************************/
void GSkyDir::lb_deg(const double& l, const double& b)
{
    // Set attributes
    m_has_lb    = 1;
    m_has_radec = 0;

    // Set direction
    m_l = l * deg2rad;
    m_b = b * deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns galactic longitude in radians
 ***************************************************************************/
double GSkyDir::l(void)
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic longitude
    return m_l;
}


/***********************************************************************//**
 * @brief Returns galactic longitude in degrees
 ***************************************************************************/
double GSkyDir::l_deg(void)
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic longitude
    return m_l * rad2deg;
}


/***********************************************************************//**
 * @brief Returns galactic latitude in radians
 ***************************************************************************/
double GSkyDir::b(void)
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic latitude
    return m_b;
}


/***********************************************************************//**
 * @brief Returns galactic latitude in degrees
 ***************************************************************************/
double GSkyDir::b_deg(void)
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic latitude
    return m_b * rad2deg;
}


/***********************************************************************//**
 * @brief Returns Right Ascension in radians
 ***************************************************************************/
double GSkyDir::ra(void)
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Right Ascension
    return m_ra;
}


/***********************************************************************//**
 * @brief Returns Right Ascension in degrees
 ***************************************************************************/
double GSkyDir::ra_deg(void)
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Right Ascension
    return m_ra * rad2deg;
}


/***********************************************************************//**
 * @brief Returns Declination in radians
 ***************************************************************************/
double GSkyDir::dec(void)
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Declination
    return m_dec;
}


/***********************************************************************//**
 * @brief Returns Declination in degrees
 ***************************************************************************/
double GSkyDir::dec_deg(void)
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Declination
    return m_dec * rad2deg;
}


/*==========================================================================
 =                                                                         =
 =                          GSkyDir private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkyDir::init_members(void)
{
    // Initialise members
    m_has_lb    = 0;
    m_has_radec = 0;
    m_l         = 0.0;
    m_b         = 0.0;
    m_ra        = 0.0;
    m_dec       = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Sky direction from which members should be copied
 ***************************************************************************/
void GSkyDir::copy_members(const GSkyDir& dir)
{
    // Copy attributes
    m_has_lb    = dir.m_has_lb;
    m_has_radec = dir.m_has_radec;
    m_l         = dir.m_l;
    m_b         = dir.m_b;
    m_ra        = dir.m_ra;
    m_dec       = dir.m_dec;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyDir::free_members(void)
{
    // Free memory

    // Signal free pointers

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert equatorial to galactic coordinates
 ***************************************************************************/
void GSkyDir::equ2gal(void)
{
    // Convert from equatorial to galactic
    euler(0, m_ra, m_dec, &m_l, &m_b);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert galactic to equatorial coordinates
 ***************************************************************************/
void GSkyDir::gal2equ(void)
{
    // Convert from galactic to equatorial
    euler(1, m_l, m_b, &m_ra, &m_dec);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GSkyDir friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] column Sky direction to put in output stream
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSkyDir& dir)
{
    // Create coordinate system dependent output
    if (dir.m_has_lb)
        os << "(l,b)=(" << dir.m_l*rad2deg << "," << dir.m_b*rad2deg << ")";
    else if (dir.m_has_radec)
        os << "(RA,Dec)=(" << dir.m_ra*rad2deg << "," << dir.m_dec*rad2deg << ")";

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GSkyDir                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief General coordinate transformation routine for J2000
 *
 * @param[in] type Conversion type (0=equ2gal, 1=gal2equ)
 * @param[in] xin Input longitude (RA or GLON) in radians.
 * @param[in] yin Input latitude (Dec or GLAT) in radians.
 * @param[out] xout Output longitude in radians.
 * @param[out] yout Output latitude in radians.
 ***************************************************************************/
void euler(const int& type, const double& xin, const double &yin, 
           double* xout, double *yout)
{
    // Set transformation constants
    const double psi[]    = {0.57477043300,  4.9368292465};
    const double stheta[] = {0.88998808748, -0.88998808748};
    const double ctheta[] = {0.45598377618,  0.45598377618};
    const double phi[]    = {4.9368292465,   0.57477043300};

    // Perform transformation
    double a    = xin - phi[type];
    double b    = yin;
    double sb   = sin(b);
    double cb   = cos(b);
    double cbsa = cb * sin(a);

    //
    a = atan2(ctheta[type] * cbsa + stheta[type] * sb, cb * cos(a));
    b = -stheta[type] * cbsa + ctheta[type] * sb;
    if (b > 1.0)
        b = 1.0;

    //
    *yout = asin(b);
    *xout = modulo((a+psi[type] + fourpi), twopi);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns the remainder of the division \a v1/v2.
 *
 * @param[in] v1 Argument 1.
 * @param[in] v2 Argument 2.
 *
 * Returns the remainder of the division \a v1/v2.
 * The result is non-negative.
 * \a v1 can be positive or negative; \a v2 must be positive.
 ***************************************************************************/
double modulo(double v1, double v2)
{
    // Return
    return (v1 >= 0) ? ((v1 < v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2);
}


/* Namespace ends ___________________________________________________________ */
}
