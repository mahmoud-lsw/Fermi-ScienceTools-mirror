/*------------------------------------------------------------------------------
Id ........: $Id: Catalogue.cxx,v 1.59 2014/06/12 12:02:35 jurgen Exp $
Author ....: $Author: jurgen $
Revision ..: $Revision: 1.59 $
Date ......: $Date: 2014/06/12 12:02:35 $
--------------------------------------------------------------------------------
$Log: Catalogue.cxx,v $
Revision 1.59  2014/06/12 12:02:35  jurgen
Added 1e-4 margin to error circle, update version number

Revision 1.58  2014/06/11 16:59:29  jurgen
Read absolute position error from 2nd FITS file extension.

Revision 1.57  2011/10/05 20:30:10  jurgen
Correctly forward row information for missing source names

Revision 1.56  2011/02/09 08:14:12  jurgen
Add "Source_Name" to possible names in counterpart catalogues (srcid.py)

Revision 1.55  2010/12/20 08:52:00  jurgen
Adapt to gcc 4.4

Revision 1.54  2010/11/04 21:20:35  jurgen
Remove unused variable

Revision 1.53  2010/10/13 14:19:35  jurgen
Dump candidates before probability cut in refine step

Revision 1.52  2010/09/26 19:01:11  jurgen
Correctly read error ellipse information.

Revision 1.51  2010/04/16 21:53:16  jurgen
Fully implement HEALPix counterpart density maps

Revision 1.50  2010/04/16 16:16:19  jurgen
Implement HEALPix interface to read counterpart density maps

Revision 1.49  2010/04/16 15:23:11  jurgen
Move euler and modula in GSkyDir.cxx

Revision 1.48  2009/03/26 14:15:05  jurgen
Properly handle NULL error radii (by assuming an minimum 1D 1sigma error of 0.005 deg)

Revision 1.47  2009/03/18 10:01:59  jurgen
Avoid floating point exception in case of NULL position errors

Revision 1.46  2008/08/20 11:52:21  jurgen
Correct probability computation and resolve STGEN-56

Revision 1.45  2008/07/08 21:21:55  jurgen
Do column name selection only from first letter on

Revision 1.44  2008/07/08 20:57:06  jurgen
Implement final selection (allows to filter on evaluated quantities)

Revision 1.43  2008/05/06 16:00:04  jurgen
Import srcid.py script and classes definition into CVS

Revision 1.42  2008/04/24 14:55:17  jurgen
Implement simple FoM scheme

Revision 1.41  2008/04/23 15:42:06  jurgen
Don't close in-memory catalogue (error if catalogue is empty)

Revision 1.40  2008/04/18 20:50:33  jurgen
Implement catch-22 scheme for prior probability calculation and compute log likelihood-ratio instead of likelihood ratio (avoid numerical problems)

Revision 1.39  2008/04/18 16:14:16  jurgen
Add LR statistics to log file

Revision 1.38  2008/04/18 10:43:20  jurgen
Allow for divergent LRs (flag them)

Revision 1.37  2008/04/16 22:00:34  jurgen
Compute unique posterior probabilities

Revision 1.36  2008/04/15 22:30:54  jurgen
Cleanup counterpart statistics

Revision 1.35  2008/04/15 21:24:12  jurgen
Introduce sparse matrix for source catalogue probability computation.

Revision 1.34  2008/04/04 14:55:52  jurgen
Remove counterpart candidate working memory and introduce permanent counterpart candidate memory

Revision 1.33  2008/03/26 16:46:58  jurgen
add more information to FITS file header

Revision 1.32  2008/03/26 13:37:10  jurgen
Generalize probability calculation and implement Bayesian method

Revision 1.31  2008/03/21 16:46:14  jurgen
Remove double log

Revision 1.30  2008/03/21 16:42:56  jurgen
Update documentation

Revision 1.29  2008/03/21 15:27:03  jurgen
Estimate number of false associations

Revision 1.28  2008/03/21 09:10:12  jurgen
Enhance code documentation.

Revision 1.27  2008/03/20 21:56:26  jurgen
implement local counterpart density

Revision 1.26  2008/03/20 12:17:44  jurgen
Invert _RA/_DE and RA/DE column name search

Revision 1.25  2008/02/23 10:36:57  jurgen
remove redundant catalogAccess header inclusion

Revision 1.24  2007/12/06 16:42:15  jurgen
Add RA/DEC and PosErr generic names

Revision 1.23  2007/11/30 16:19:26  jurgen
Correct version number and add RAdeg/DEdeg columns

Revision 1.22  2007/11/08 14:42:11  jurgen
Handle error circles (e.g. 3EG catalogue)

Revision 1.21  2007/10/11 13:20:54  jurgen
Correctly remove FITS special function columns

Revision 1.20  2007/10/10 15:39:12  jurgen
Introduce handling of special functions 'gammln', 'erf', and 'erfc'

Revision 1.19  2007/10/09 16:46:23  jurgen
Write counterpart catalogue reference (row) to output catalogue

Revision 1.18  2007/10/09 08:17:40  jurgen
Correctly interpret positional errors and correctly evaluate PROB_POS
as likelihood

Revision 1.17  2007/10/08 11:02:25  jurgen
Implement search for catalogue table information and handle different
position error types

Revision 1.16  2007/09/21 20:27:14  jurgen
Correct cfits_collect bug (unstable row selection)

Revision 1.15  2007/09/21 14:29:03  jurgen
Correct memory bug and updated test script

Revision 1.14  2007/09/21 12:49:10  jurgen
Enhance log-file output and chatter level

Revision 1.13  2007/09/20 16:28:21  jurgen
Enhance catalogue interface for column recognition

Revision 1.12  2006/02/09 13:06:18  jurgen
Put maximum number of source to load at once in a constant and change
value to a large value (since the loading logic has not yet been
implemented).

Revision 1.11  2006/02/07 16:05:04  jurgen
Use ObjectInfo structure to hold catalogue object information

Revision 1.10  2006/02/07 11:10:50  jurgen
Suppress catalogAccess verbosity

Revision 1.9  2006/02/03 22:10:30  jurgen
Remove comments to correctly catch an error in the determination
of the number of catalogue entries (these comments have been
introduced due to a bug in catalogAccess - which has been fixed in
version v0r2p6).

Revision 1.8  2006/02/03 12:14:51  jurgen
New version that allows additional probabilities to be taken
into account. The code has been considerably reorganised. Also
catalogue column prefixes are now handled differently.

Revision 1.7  2006/02/02 09:31:55  jurgen
fix last Win32 bug !!!

Revision 1.6  2006/02/02 09:29:29  jurgen
correct Win32 bug

Revision 1.5  2006/02/02 09:26:14  jurgen
correct Win32 compile bugs

Revision 1.4  2006/02/01 15:17:10  jurgen
correct g++-3.4.3 compile error related to parentheses around (char*)

Revision 1.3  2006/02/01 13:33:36  jurgen
Tried to fix Win32 compilation bugs.
Change revision number to 1.3.2.
Replace header information with CVS typeset information.

------------------------------------------------------------------------------*/
/**
 * @file Catalogue.cxx
 * @brief Implements methods of Catalogue class.
 * @author J. Knodlseder
 */

/* Includes _________________________________________________________________ */
#include <exception>
#include <cstring>
#include "sourceIdentify.h"
#include "Catalogue.h"
#include "Log.h"


/* Definitions ______________________________________________________________ */
#define CATALOGUE_TIMING   0                     // Enables timing measurements


/* Namespace definition _____________________________________________________ */
namespace sourceIdentify {
using namespace catalogAccess;


/* Globals __________________________________________________________________ */


/* Constants ________________________________________________________________ */


/* Type defintions __________________________________________________________ */


/* Private Prototypes _______________________________________________________ */
std::string find(std::vector <std::string> &arg, std::string match);
double      modulo(double v1, double v2);
void        euler(const int& type, const double& xin, const double &yin, 
                  double* xout, double *yout);
Status      get_info(Parameters *par, InCatalogue *in, Status status);
Status      get_id_info(Parameters *par, InCatalogue *in, 
                        std::vector <std::string> &qtyNames,
                        std::vector <std::string> &qtyUCDs,
                        Status status);
Status      get_pos_info(Parameters *par, InCatalogue *in,
                         std::vector <std::string> &qtyNames,
                         std::vector <std::string> &qtyUCDs,
                         Status status);
Status      get_pos_error_info(Parameters *par, InCatalogue *in,
                               std::vector <std::string> &qtyNames,
                               Status status);
void        set_info(Parameters *par, InCatalogue *in, int &i, ObjectInfo *ptr,
                     double &posErr);


/*============================================================================*/
/*                              Private functions                             */
/*============================================================================*/

/**************************************************************************//**
 * @brief Convert string to upper case
 *
 * @param[in] arg String to convert into upper case.
 ******************************************************************************/
std::string upper(std::string arg) {

  // Copy argument
  std::string result = arg;

  // Convert to upper case
  std::transform(result.begin(), result.end(), result.begin(),
                 (int(*)(int)) std::toupper);

  // Return result
  return result;

}


/**************************************************************************//**
 * @brief Returns the shortest string that matches string
 *
 * @param[in] arg Vector of strings to be compated to match.
 * @param[in] match String to be matched.
 *
 * From a list of strings, returns the string that matches the pattern specified
 * by 'match'. If more than a single match exists, the shortest of all matches
 * is returned.
 ******************************************************************************/
std::string find(std::vector <std::string> &arg, std::string match) {

  // Initialise result
  std::string result;
  int         length = 0;

  // Convert match to upper case
  match = upper(match);

  // Loop over all vector elements
  for (int i = 0; i < (int)arg.size(); ++i) {
//    if (upper(arg[i]).find(match, 0) != std::string::npos) {
    if (upper(arg[i]).find(match, 0) == 0) {
      int this_length = arg[i].length();
      if (length == 0) {
        length = this_length;
        result = arg[i];
      }
      else if (this_length < length) {
        length = this_length;
        result = arg[i];
      }
    }
  }

  // Make sure that result string is empty if we found nothing
  if (length == 0)
    result.clear();

  // Return result
  return result;

}


/**************************************************************************//**
 * @brief Get column information
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] in Pointer to input catalogue.
 * @param[in] status Error status.
 ******************************************************************************/
Status get_info(Parameters *par, InCatalogue *in, Status status) {

    // Declare local variables
    int                       numKeys;
    int                       numUCDs;
    std::vector <std::string> qtyNames;
    std::vector <std::string> qtyUCDs;

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: get_info");

    // Single loop for common exit point
    do {

      // Fall through in case of an error
      if (status != STATUS_OK)
        continue;

      // Extract keys and UCDs from catalogue. Stop if their number does not
      // correspond
      numKeys = in->cat.getQuantityNames(&qtyNames);
      numUCDs = in->cat.getQuantityUCDs(&qtyUCDs);
      if (numKeys != numUCDs) {
        if (par->logTerse())
          Log(Error_2, "%d : Mismatch between number of keys (%d) and UCDs (%d)",
              status, numKeys, numUCDs);
        continue;
      }

      // Get source ID column
      status = get_id_info(par, in, qtyNames, qtyUCDs, status);
      if (status == STATUS_CAT_NO_ID) {
        status = STATUS_OK;
        in->col_id.clear();
      }

      // Get position columns
      status = get_pos_info(par, in, qtyNames, qtyUCDs, status);
      if (status == STATUS_CAT_NO_POS) {
        status = STATUS_OK;
        in->col_ra.clear();
        in->col_dec.clear();
      }

      // Get position error columns
      status = get_pos_error_info(par, in, qtyNames, status);
      if (status == STATUS_CAT_NO_POS_ERROR) {
        status         = STATUS_OK;
        in->col_e_type = NoError;
      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: get_info (status=%d)", 
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Get source name information
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] in Pointer to input catalogue.
 * @param[in] status Error status.
 *
 * The following ranked criteria are applied to determine the name of the
 * catalogue column that contains the source name:
 * 1) Search for column with UCD 'ID_MAIN'
 * 2) Search for column with names contained in search string (NAME, ID)
 ******************************************************************************/
Status get_id_info(Parameters *par, InCatalogue *in,
                   std::vector <std::string> &qtyNames,
                   std::vector <std::string> &qtyUCDs,
                   Status status) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: get_id_info");

    // Single loop for common exit point
    do {

      // Fall through in case of an error
      if (status != STATUS_OK)
        continue;

      // Initialise status to 'not found'
      status = STATUS_CAT_NO_ID;

      // Search for first ID in UCDs
      for (int i = 0; i < (int)qtyUCDs.size(); ++i) {
        if (qtyUCDs[i].find("ID_MAIN", 0) != std::string::npos) {
          in->col_id = qtyNames[i];
          status     = STATUS_OK;
          break;
        }
      }
      if (status == STATUS_OK)
        continue;

      // Search for ID in keys
      for (int k = 0; search_id[k] != "stop"; ++k) {
        std::string match = find(qtyNames, upper(search_id[k]));
        if (match.length() > 0) {
          in->col_id = match;
          status     = STATUS_OK;
          break;
        }
      }
      if (status == STATUS_OK)
        continue;

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: get_id_info (status=%d)", 
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Get source position information
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] in Pointer to input catalogue.
 * @param[in] status Error status.
 *
 * The following ranked criteria are applied to determine the name of the
 * catalogue columns that contain the source position:
 * 1) Search for columns with UCDs 'POS_EQ_RA_MAIN' & 'POS_EQ_DEC_MAIN'
 * 2) Search for columns with names 'RAdeg' & 'DEdeg'
 * 3) Search for columns with names '_RAJ2000' & '_DEJ2000'
 * 4) Search for columns with names 'RAJ2000' & 'DEJ2000'
 * 5) Search for columns with names 'RA' & 'DEC'
 ******************************************************************************/
Status get_pos_info(Parameters *par, InCatalogue *in,
                    std::vector <std::string> &qtyNames,
                    std::vector <std::string> &qtyUCDs,
                    Status status) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: get_pos_info");

    // Single loop for common exit point
    do {

      // Fall through in case of an error
      if (status != STATUS_OK)
        continue;

      // Initialise status to 'not found'
      status       = STATUS_CAT_NO_POS;
      in->pos_type = NoPosition;
      in->col_ra.clear();
      in->col_dec.clear();
      in->col_glon.clear();
      in->col_glat.clear();

      // Search for first RA/Dec position in UCDs
      for (int i = 0; i < (int)qtyUCDs.size(); ++i) {
        if (qtyUCDs[i].find("POS_EQ_RA_MAIN", 0) != std::string::npos) {
          in->col_ra   = qtyNames[i];
          break;
        }
      }
      for (int i = 0; i < (int)qtyUCDs.size(); ++i) {
        if (qtyUCDs[i].find("POS_EQ_DEC_MAIN", 0) != std::string::npos) {
          in->col_dec  = qtyNames[i];
          break;
        }
      }
      if ((in->col_ra.length() > 0) && (in->col_dec.length() > 0)) {
        status       = STATUS_OK;
        in->pos_type = Equatorial;
        continue;
      }

      // Declare column string
      std::string col_ra;
      std::string col_dec;
      std::string col_glon;
      std::string col_glat;

      // Search for RAdeg/DEdeg columns
      col_ra  = find(qtyNames, "RAdeg");
      col_dec = find(qtyNames, "DEdeg");
      if ((col_ra.length() > 0) && (col_dec.length() > 0)) {
        in->col_ra   = col_ra;
        in->col_dec  = col_dec;
        in->pos_type = Equatorial;
        status       = STATUS_OK;
        continue;
      }

      // Search for _RAJ2000/_DEJ2000 columns
      col_ra  = find(qtyNames, "_RAJ2000");
      col_dec = find(qtyNames, "_DEJ2000");
      if ((col_ra.length() > 0) && (col_dec.length() > 0)) {
        in->col_ra   = col_ra;
        in->col_dec  = col_dec;
        in->pos_type = Equatorial;
        status       = STATUS_OK;
        continue;
      }

      // Search for RAJ2000/DEJ2000 columns
      col_ra  = find(qtyNames, "RAJ2000");
      col_dec = find(qtyNames, "DEJ2000");
      if ((col_ra.length() > 0) && (col_dec.length() > 0)) {
        in->col_ra   = col_ra;
        in->col_dec  = col_dec;
        in->pos_type = Equatorial;
        status       = STATUS_OK;
        continue;
      }

      // Search for RA/DEC columns
      col_ra  = find(qtyNames, "RA");
      col_dec = find(qtyNames, "DEC");
      if ((col_ra.length() > 0) && (col_dec.length() > 0)) {
        in->col_ra   = col_ra;
        in->col_dec  = col_dec;
        in->pos_type = Equatorial;
        status       = STATUS_OK;
        continue;
      }

      // Search for GLON/GLAT position in UCDs
      for (int i = 0; i < (int)qtyUCDs.size(); ++i) {
        if (qtyUCDs[i].find("POS_GAL_LON", 0) != std::string::npos) {
          in->col_glon = qtyNames[i];
          break;
        }
      }
      for (int i = 0; i < (int)qtyUCDs.size(); ++i) {
        if (qtyUCDs[i].find("POS_GAL_LAT", 0) != std::string::npos) {
          in->col_glat = qtyNames[i];
          break;
        }
      }
      if ((in->col_glon.length() > 0) && (in->col_glat.length() > 0)) {
        status       = STATUS_OK;
        in->pos_type = Galactic;
        continue;
      }

      // Search for _GLON/_GLAT columns
      col_glon = find(qtyNames, "_GLON");
      col_glat = find(qtyNames, "_GLAT");
      if ((col_glon.length() > 0) && (col_glat.length() > 0)) {
        in->col_glon = col_glon;
        in->col_glat = col_glat;
        in->pos_type = Galactic;
        status       = STATUS_OK;
        continue;
      }

      // Search for GLON/GLAT columns
      col_glon = find(qtyNames, "GLON");
      col_glat = find(qtyNames, "GLAT");
      if ((col_glon.length() > 0) && (col_glat.length() > 0)) {
        in->col_glon = col_glon;
        in->col_glat = col_glat;
        in->pos_type = Galactic;
        status       = STATUS_OK;
        continue;
      }

      // Search for L/B columns
      col_glon = find(qtyNames, "L");
      col_glat = find(qtyNames, "B");
      if ((col_glon.length() > 0) && (col_glat.length() > 0)) {
        in->col_glon = col_glon;
        in->col_glat = col_glat;
        in->pos_type = Galactic;
        status       = STATUS_OK;
        continue;
      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: get_pos_info (status=%d)", status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Get source position error information
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] in Pointer to input catalogue.
 * @param[in] status Error status.
 *
 * The following ranked criteria are applied to determine the name of the
 * catalogue columns that contain the source position:
 * 1) Search for columns with names 'Conf_95_SemiMajor', 'Conf_95_SemiMinor' &
 *    'Conf_95_PosAng'
 * 2) Search for columns with names 'Conf_68_SemiMajor', 'Conf_68_SemiMinor' &
 *    'Conf_68_PosAng'
 * 3) Search for columns with names OUTCAT_COL_MAJERR_NAME, 
 *    OUTCAT_COL_MINERR_NAME & OUTCAT_COL_POSANGLE_NAME
 * 4) Search for columns with names 'e_RAdeg' & 'e_DEdeg'
 * 5) Search for columns with names 'e_RAJ2000' & 'e_DEJ2000'
 * 6) Search for columns with name 'theta95'
 * 7) Search for columns with name 'PosErr'
 ******************************************************************************/
Status get_pos_error_info(Parameters *par, InCatalogue *in,
                          std::vector <std::string> &qtyNames,
                          Status status) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: get_pos_error_info");

    // Single loop for common exit point
    do {

      // Fall through in case of an error
      if (status != STATUS_OK)
        continue;

      // Initialise status to 'not found'
      status = STATUS_CAT_NO_POS_ERROR;

      // Initialise error type to 'no error' and error scaling to unity
      in->col_e_type  = NoError;
      in->e_pos_scale = 1.0;

      // Search loop
      do {

        // Declare column string
        std::string col_e_maj;
        std::string col_e_min;
        std::string col_e_posang;
        std::string col_e_ra;
        std::string col_e_dec;

        // Search for LAT catalogue names (95%)
        col_e_maj    = find(qtyNames, "Conf_95_SemiMajor");
        col_e_min    = find(qtyNames, "Conf_95_SemiMinor");
        col_e_posang = find(qtyNames, "Conf_95_PosAng");
        if ((col_e_maj.length() > 0) && (col_e_min.length() > 0) &&
            (col_e_posang.length() > 0)) {
          in->col_e_maj    = col_e_maj;
          in->col_e_min    = col_e_min;
          in->col_e_posang = col_e_posang;
          in->col_e_type   = Ellipse;
          in->col_e_prob   = Prob_95;
          status           = STATUS_OK;
          continue;
        }

        // Search for LAT catalogue names (68%)
        col_e_maj    = find(qtyNames, "Conf_68_SemiMajor");
        col_e_min    = find(qtyNames, "Conf_68_SemiMinor");
        col_e_posang = find(qtyNames, "Conf_68_PosAng");
        if ((col_e_maj.length() > 0) && (col_e_min.length() > 0) &&
            (col_e_posang.length() > 0)) {
          in->col_e_maj    = col_e_maj;
          in->col_e_min    = col_e_min;
          in->col_e_posang = col_e_posang;
          in->col_e_type   = Ellipse;
          in->col_e_prob   = Prob_68;
          status           = STATUS_OK;
          continue;
        }

        // Search for output catalogue columns
        col_e_maj    = find(qtyNames, OUTCAT_COL_MAJERR_NAME);
        col_e_min    = find(qtyNames, OUTCAT_COL_MINERR_NAME);
        col_e_posang = find(qtyNames, OUTCAT_COL_POSANGLE_NAME);
        if ((col_e_maj.length() > 0) && (col_e_min.length() > 0) &&
            (col_e_posang.length() > 0)) {
          in->col_e_maj    = col_e_maj;
          in->col_e_min    = col_e_min;
          in->col_e_posang = col_e_posang;
          in->col_e_type   = Ellipse;
          in->col_e_prob   = Prob_95;
          status           = STATUS_OK;
          continue;
        }

        // Search for e_RAdeg/e_DEdeg columns
        if (in->pos_type == Equatorial) {
          col_e_ra  = find(qtyNames, "e_RAdeg");
          col_e_dec = find(qtyNames, "e_DEdeg");
          if ((col_e_ra.length() > 0) && (col_e_dec.length() > 0)) {
            in->col_e_ra   = col_e_ra;
            in->col_e_dec  = col_e_dec;
            in->col_e_type = RaDec;
            in->col_e_prob = Sigma_1;
            status         = STATUS_OK;
            continue;
          }
        }

        // Search for e_RAJ2000/e_DEJ2000 columns
        if (in->pos_type == Equatorial) {
          col_e_ra  = find(qtyNames, "e_RAJ2000");
          col_e_dec = find(qtyNames, "e_DEJ2000");
          if ((col_e_ra.length() > 0) && (col_e_dec.length() > 0)) {
            in->col_e_ra   = col_e_ra;
            in->col_e_dec  = col_e_dec;
            in->col_e_type = RaDec;
            in->col_e_prob = Sigma_1;
            status         = STATUS_OK;
            continue;
          }
        }

        // Search for theta95 column (3EG catalogue)
        col_e_maj = find(qtyNames, "theta95");
        if (col_e_maj.length() > 0) {
          in->col_e_maj  = col_e_maj;
          in->col_e_type = Radius;
          in->col_e_prob = Prob_95;
          status         = STATUS_OK;
          continue;
        }

        // Search for PosErr68 column
        col_e_maj = find(qtyNames, "PosErr68");
        if (col_e_maj.length() > 0) {
          in->col_e_maj  = col_e_maj;
          in->col_e_type = Radius;
          in->col_e_prob = Prob_68;
          status         = STATUS_OK;
          continue;
        }

        // Search for PosErr90 column
        col_e_maj = find(qtyNames, "PosErr90");
        if (col_e_maj.length() > 0) {
          in->col_e_maj  = col_e_maj;
          in->col_e_type = Radius;
          in->col_e_prob = Prob_90;
          status         = STATUS_OK;
          continue;
        }

        // Search for PosErr95 column
        col_e_maj = find(qtyNames, "PosErr95");
        if (col_e_maj.length() > 0) {
          in->col_e_maj  = col_e_maj;
          in->col_e_type = Radius;
          in->col_e_prob = Prob_95;
          status         = STATUS_OK;
          continue;
        }

        // Search for PosErr99 column
        col_e_maj = find(qtyNames, "PosErr99");
        if (col_e_maj.length() > 0) {
          in->col_e_maj  = col_e_maj;
          in->col_e_type = Radius;
          in->col_e_prob = Prob_99;
          status         = STATUS_OK;
          continue;
        }

        // Search for PosErr column
        col_e_maj = find(qtyNames, "PosErr");
        if (col_e_maj.length() > 0) {
          in->col_e_maj  = col_e_maj;
          in->col_e_type = Radius;
          in->col_e_prob = Sigma_1;
          status         = STATUS_OK;
          continue;
        }

      } while (0); // End of search loop

      // Get error scaling (returned errors are 95% confidence errors)
      if (in->col_e_type != NoError) {
        switch (in->col_e_prob) {
        case Sigma_1:
          in->e_pos_scale = e_norm_1s;
          break;
        case Sigma_2:
          in->e_pos_scale = e_norm_2s;
          break;
        case Sigma_3:
          in->e_pos_scale = e_norm_3s;
          break;
        case Prob_68:
          in->e_pos_scale = e_norm_68;
          break;
        case Prob_90:
          in->e_pos_scale = e_norm_90;
          break;
        case Prob_95:
          in->e_pos_scale = e_norm_95;
          break;
        case Prob_99:
          in->e_pos_scale = e_norm_99;
          break;
        default:
          in->e_pos_scale = e_norm_95;
          break;
        }
        in->e_pos_scale /= e_norm_95;
      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: get_pos_error_info (status=%d)", status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Set information for source from catalogue
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] in Pointer to input catalogue.
 * @param[in] i Source number (starting from 0).
 * @param[in] ptr Pointer to source information structure.
 * @param[in] posErr Error radius if no error is found in catalogue.
 ******************************************************************************/
void set_info(Parameters *par, InCatalogue *in, int &i, ObjectInfo *ptr,
              double &posErr) {

    // Declare local variables
    double err_maj;
    double err_min;
    double err_ang;
    double e_RA;
    double e_DE;

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: set_info");

    // Single loop for common exit point
    do {

      // Initialise source information
      ptr->pos_valid   = 0;        // Invalid position
      ptr->name.clear();
      ptr->pos_eq_ra   = 0.0;
      ptr->pos_eq_dec  = 0.0;
      ptr->pos_err_maj = posErr;
      ptr->pos_err_min = posErr;
      ptr->pos_err_ang = 0.0;

      // Set source name
      if (in->cat.getSValue(in->col_id, i, &(ptr->name)) != IS_OK)
        ptr->name = "no-name";

      // Set source position
      if (in->pos_type == Equatorial) {
        if ((in->cat.getNValue(in->col_ra,  i, &(ptr->pos_eq_ra))  == IS_OK) &&
            (in->cat.getNValue(in->col_dec, i, &(ptr->pos_eq_dec)) == IS_OK)) {

          // Set position validity flag
          ptr->pos_valid = 1;

          // Put Right Ascension in interval [0,2pi[
          ptr->pos_eq_ra = ptr->pos_eq_ra -
                           double(long(ptr->pos_eq_ra / 360.0) * 360.0);
          if (ptr->pos_eq_ra < 0.0)
            ptr->pos_eq_ra += 360.0;

        } // endif: equatorial source position found
      }
      else if (in->pos_type == Galactic) {
        double glon;
        double glat;
        if ((in->cat.getNValue(in->col_glon, i, &glon)  == IS_OK) &&
            (in->cat.getNValue(in->col_glat, i, &glat) == IS_OK)) {

          // Convert galactic to equatorial coordinates
          glon *= deg2rad;
          glat *= deg2rad;
          euler(1, glon, glat, &ptr->pos_eq_ra, &ptr->pos_eq_dec);
          ptr->pos_eq_ra  *= rad2deg;
          ptr->pos_eq_dec *= rad2deg;

          // Set position validity flag
          ptr->pos_valid = 1;

          // Put Right Ascension in interval [0,2pi[
          ptr->pos_eq_ra = ptr->pos_eq_ra -
                           double(long(ptr->pos_eq_ra / 360.0) * 360.0);
          if (ptr->pos_eq_ra < 0.0)
            ptr->pos_eq_ra += 360.0;

        } // endif: galactic source position found
      }

      // Set source position error (type dependent)
      switch (in->col_e_type) {
      case NoError:
        ptr->pos_err_maj = posErr;
        ptr->pos_err_min = posErr;
        ptr->pos_err_ang = 0.0;
        break;
      case Radius:
        if (in->cat.getNValue(in->col_e_maj, i, &err_maj) == IS_OK) {
          ptr->pos_err_maj = err_maj * in->e_pos_scale;
          ptr->pos_err_min = err_maj * in->e_pos_scale;
          ptr->pos_err_ang = 0.0;
        }
        break;
      case Ellipse:
        if ((in->cat.getNValue(in->col_e_maj,    i, &err_maj) == IS_OK) &&
            (in->cat.getNValue(in->col_e_min,    i, &err_min) == IS_OK) &&
            (in->cat.getNValue(in->col_e_posang, i, &err_ang) == IS_OK)) {
          ptr->pos_err_maj = err_maj * in->e_pos_scale;
          ptr->pos_err_min = err_min * in->e_pos_scale;
          ptr->pos_err_ang = err_ang;
        }
        break;
      case RaDec:
        if ((in->cat.getNValue(in->col_e_ra, i, &e_RA) == IS_OK) &&
            (in->cat.getNValue(in->col_e_dec, i, &e_DE) == IS_OK)) {
          e_RA *= cos(ptr->pos_eq_dec*deg2rad);
          if (e_RA > e_DE) {           // Error ellipse along RA axis
            ptr->pos_err_maj = e_RA * in->e_pos_scale;
            ptr->pos_err_min = e_DE * in->e_pos_scale;
            ptr->pos_err_ang = 90.0;   // P.A. = 90.0 deg
          }
          else {                       // Error ellipse along DE axis
            ptr->pos_err_maj = e_DE * in->e_pos_scale;
            ptr->pos_err_min = e_RA * in->e_pos_scale;
            ptr->pos_err_ang = 0.0;    // P.A. = 0.0 deg
          }
        }
        break;
      }

      // Avoid source positions errors smaller than the absolute position error
      if (ptr->pos_err_maj < in->erposabs && ptr->pos_err_min < in->erposabs) {
        ptr->pos_err_maj = in->erposabs;
        ptr->pos_err_min = in->erposabs;
        ptr->pos_err_ang = 0.0;
      }

     } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: set_info");

    // Return
    return;

}


/*============================================================================*/
/*                          Low-level catalogue methods                       */
/*============================================================================*/

/**************************************************************************//**
 * @brief Initialise class memory.
 *
 * The maximum number of counterparts is set by the constant c_maxCptLoad.
 * The filter step bounding box size is set by the constant c_filter_maxsep.
 ******************************************************************************/
void Catalogue::init_memory(void) {

    // Declare local variables

    // Single loop for common exit point
    do {

      // Initialise source catalogue private members
      m_src.numLoad     = 0;
      m_src.numTotal    = 0;
      m_src.object      = NULL;
      m_src.col_e_type  = NoError;
      m_src.e_pos_scale = 1.0;
      m_src.inName.clear();
      m_src.catCode.clear();
      m_src.catURL.clear();
      m_src.catName.clear();
      m_src.catRef.clear();
      m_src.tableName.clear();
      m_src.tableRef.clear();
      m_src.col_id.clear();
      m_src.col_ra.clear();
      m_src.col_dec.clear();
      m_src.col_e_ra.clear();
      m_src.col_e_dec.clear();
      m_src.col_e_maj.clear();
      m_src.col_e_min.clear();

      // Initialise counterpart catalogue private members
      m_cpt.numLoad     = 0;
      m_cpt.numTotal    = 0;
      m_cpt.object      = NULL;
      m_cpt.col_e_type  = NoError;
      m_cpt.e_pos_scale = 1.0;
      m_cpt.inName.clear();
      m_cpt.catCode.clear();
      m_cpt.catURL.clear();
      m_cpt.catName.clear();
      m_cpt.catRef.clear();
      m_cpt.tableName.clear();
      m_cpt.tableRef.clear();
      m_cpt.col_id.clear();
      m_cpt.col_ra.clear();
      m_cpt.col_dec.clear();
      m_cpt.col_e_ra.clear();
      m_cpt.col_e_dec.clear();
      m_cpt.col_e_maj.clear();
      m_cpt.col_e_min.clear();

      // Intialise catalogue FITS files
      m_memFile = NULL;
      m_outFile = NULL;

      // Initialise source information
      m_info = NULL;

      // Initialise list of selected counterparts
      m_cpt_sel = NULL;

      // Initialise counterpart statistics
      m_num_Sel  = 0;
      m_cpt_stat = NULL;

      // Initialise counterpart density flag
      m_has_density = 0;

      // Catch-22
      m_prior     = c_prob_prior;
      m_prior_min = c_prob_prior_min;
      m_prior_max = c_prob_prior_max;
      m_iter      = 0;

      // Initialise association results
      m_num_claimed      = 0.0;
      m_sum_pid          = 0.0;
      m_sum_pc           = 0.0;
      m_sum_lr           = 0.0;
      m_sum_pid_thr      = 0.0;
      m_sum_pc_thr       = 0.0;
      m_sum_lr_thr       = 0.0;
      m_reliability      = 0.0;
      m_completeness     = 0.0;
      m_fract_not_unique = 0.0;
      m_num_lr_div       = 0.0;

      // Initialise output catalogue quantities
      m_num_src_Qty   = 0;
      m_num_cpt_Qty   = 0;

    } while (0); // End of main do-loop

    // Return
    return;

}


/**************************************************************************//**
 * @brief Free class memory.
 ******************************************************************************/
void Catalogue::free_memory(void) {

    // Declare local variables

    // Single loop for common exit point
    do {

      // Free source information
      if (m_info != NULL) {
        for (int i = 0; i < m_src.numLoad; ++i) {
          if (m_info[i].cc != NULL) delete [] m_info[i].cc;
        }
        delete [] m_info;
      }

      // Free temporary memory
      if (m_src.object != NULL) delete [] m_src.object;
      if (m_cpt.object != NULL) delete [] m_cpt.object;
      if (m_cpt_stat   != NULL) delete [] m_cpt_stat;
      if (m_cpt_sel    != NULL) delete [] m_cpt_sel;

      // Initialise memory
      init_memory();

    } while (0); // End of main do-loop

    // Return
    return;

}


/*============================================================================*/
/*                         High-level catalogue methods                       */
/*============================================================================*/

/**************************************************************************//**
 * @brief Get descriptor for input catalogue
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] catName Catalogue name.
 * @param[in] in Pointer to input catalogue.
 * @param[in] status Error status.
 ******************************************************************************/
Status Catalogue::get_input_descriptor(Parameters *par, std::string catName,
                                       InCatalogue *in, Status status) {

    // Declare local variables
    int                      caterr;
    std::vector<std::string> titles;
    int                      fstatus;
    //char                     comment[80];
    fitsfile                *fptr;


    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::get_input_descriptor");

    // Single loop for common exit point
    do {

      // Set catalogAccess verbosity
      catalogAccess::verbosity = g_u9_verbosity;

      // Determine the number of objects in the catalogue. First we try to
      // access the catalogue on disk, then on the Web ...
      caterr = in->cat.getMaxNumRows(&in->numTotal, catName);
      if (caterr < 0) {
        if (par->logVerbose())
          Log(Warning_2, "%d : Unable to determine catalogue '%s' size from"
              " file. Try on Web now.", caterr, catName.c_str());
        caterr = in->cat.getMaxNumRowsWeb(&in->numTotal, catName);
        if (caterr < 0) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to determine number of objects in"
                " catalogue '%s'.", caterr, catName.c_str());
          status = STATUS_CAT_NOT_FOUND;
          continue;
        }
      }

      // Import the catalogue descriptor. First we try to access the catalogue 
      // on disk, then on the Web ...
      caterr = in->cat.importDescription(catName);
      if (caterr < 0) {
        if (par->logVerbose())
          Log(Warning_2, "%d : Unable to load catalogue '%s' descriptor from"
              " file. Try on Web now.", caterr, catName.c_str());
        caterr = in->cat.importDescriptionWeb(catName);
        if (caterr < 0) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to load catalogue '%s' descriptor from"
                " file or web.", caterr, catName.c_str());
          status = STATUS_CAT_NOT_FOUND;
          continue;
        }
        else {
          // Acknowledge loading
          if (par->logVerbose())
            Log(Log_2, " Loaded catalogue '%s' descriptor from web.",
                catName.c_str());

          // Set ERPOSABS keyword in 2D 95% units
          in->erposabs = c_erposabs * 2.4860;
        }
      }
      else {
        // Acknowledge loading
        if (par->logVerbose())
          Log(Log_2, " Loaded catalogue '%s' descriptor from file.",
              catName.c_str());

        // Try reading ERPOSABS keyword from 2nd extension and convert in 2D 95% units
        /*
        fstatus = 0;
        fstatus = fits_open_file(&fptr, catName.c_str(), 0, &fstatus);
        fstatus = fits_movabs_hdu(fptr, 2, NULL, &fstatus);
        fstatus = fits_read_key_dbl(fptr, "ERPOSABS", &in->erposabs, NULL, &fstatus);
        fstatus = fits_close_file(fptr, &fstatus);
        if (fstatus != 0) {
          in->erposabs = c_erposabs;
          if (par->logTerse()) {
            Log(Warning_2, " No ERPOSABS keyword found in catalogue %s (status=%d).",
                catName.c_str(), fstatus);
            Log(Warning_2, " Assume ERPOSABS value of %f.", in->erposabs);
          }
        }
        */

        // Set ERPOSABS keyword in 2D 95% units
        // Do no longer read the keywords as the catalogue applies the corrections directly
        in->erposabs = c_erposabs * 2.4860;
      }

      // Store input name
      in->inName = catName;

      // Set title vector to 6 elements
      titles.clear();
      titles.push_back(" ");
      titles.push_back(" ");
      titles.push_back(" ");
      titles.push_back(" ");
      titles.push_back(" ");
      titles.push_back(" ");

      // Extract titles from catalogue
      in->cat.getCatalogTitles(&titles);

      // Store catalogue title information
      in->catCode   = titles[0];
      in->catURL    = titles[1];
      in->catName   = titles[2];
      in->catRef    = titles[3];
      in->tableName = titles[4];
      in->tableRef  = titles[5];

      // Determine the number of loaded objects in catalogue (should be 0)
      in->cat.getNumRows(&in->numLoad);

      // Get information
      status = get_info(par, in, status);
      if (status != STATUS_OK)
        continue;

      // Dump catalogue descriptor (optionally)
//      if (par->logTerse()) {
//        status = dump_descriptor(par, in, status);
//        if (status != STATUS_OK)
//          continue;
//      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::get_input_descriptor (status=%d)", 
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Get data for input catalogue
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] in Pointer to input catalogue.
 * @param[in] posErr Error radius in case that information is missing in catalogue.
 * @param[in] status Error status.
 ******************************************************************************/
Status Catalogue::get_input_catalogue(Parameters *par, InCatalogue *in,
                                      double posErr, Status status) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::get_input_catalogue");

    // Single loop for common exit point
    do {

      // Set catalogAccess verbosity
      catalogAccess::verbosity = g_u9_verbosity;

      // First interpret the input string as filename and load the catalogue
      // from the file. If this fails then interpret input string as
      // catalogue name and load from WEB.
      int caterr = in->cat.import(in->inName);
      if (caterr < 0) {
        if (par->logVerbose())
          Log(Warning_2, "%d : Unable to load catalogue '%s' from file.",
              caterr, in->inName.c_str());
        caterr = in->cat.importWeb(in->inName);
        if (caterr < 0) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to load catalogue '%s' from file or web.",
                caterr, in->inName.c_str());
          status = STATUS_CAT_NOT_FOUND;
          continue;
        }
        else {
          if (par->logVerbose())
            Log(Log_2, " Loaded catalogue '%s' from web.", in->inName.c_str());
        }
      }
      else {
        if (par->logVerbose())
          Log(Log_2, " Loaded catalogue '%s' from file.", in->inName.c_str());
      }

      // Determine the number of loaded objects in catalogue. Fall throgh if
      // there are no objects loaded
      in->cat.getNumRows(&in->numLoad);
      if (in->numLoad < 1)
        continue;

      // Allocate memory for object information
      if (in->object != NULL) delete [] in->object;
      in->object = new ObjectInfo[in->numLoad];
      if (in->object == NULL) {
        status = STATUS_MEM_ALLOC;
        if (par->logTerse())
          Log(Error_2, "%d : Memory allocation failure.", (Status)status);
        continue;
      }

      // Extract object information
      ObjectInfo *ptr = in->object;
      for (int i = 0; i < in->numLoad; i++, ptr++) {

        // Set source information
        set_info(par, in, i, ptr, posErr);

        // Assign source name
        ptr->name = cid_assign_src_name(ptr->name, i);

      } // endfor: looped over all objects

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::get_input_catalogue (status=%d)",
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Dump catalogue descriptor
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] in Pointer to input catalogue.
 * @param[in] status Error status.
 ******************************************************************************/
Status Catalogue::dump_descriptor(Parameters *par, InCatalogue *in,
                                  Status status) {

    // Declare local variables
    long                                               iQty;
    long                                               numQty;
    long                                               len;
    long                                               maxLenNames;
    long                                               maxLenUnits;
    long                                               maxLenForms;
    long                                               maxLenUCDs;
    std::string                                        qtyType;
    std::vector<std::string>                           qtyNames;
    std::vector<std::string>                           qtyUnits;
    std::vector<std::string>                           qtyUCDs;
    std::vector<catalogAccess::Quantity>               qtyDesc;
    std::vector<catalogAccess::Quantity::QuantityType> qtyTypes;

    // Single loop for common exit point
    do {

      // Extract information from catalogue
      numQty = in->cat.getQuantityNames(&qtyNames);
      numQty = in->cat.getQuantityUnits(&qtyUnits);
      numQty = in->cat.getQuantityUCDs(&qtyUCDs);
      numQty = in->cat.getQuantityDescription(&qtyDesc);
      numQty = in->cat.getQuantityTypes(&qtyTypes);

      // Get string lengths
      maxLenNames = 0;
      maxLenUnits = 0;
      maxLenForms = 0;
      maxLenUCDs  = 0;
      for (iQty = 0; iQty < numQty; iQty++) {
        if ((len = strlen(qtyNames[iQty].c_str())) > maxLenNames)
          maxLenNames = len;
        if ((len = strlen(qtyUnits[iQty].c_str())) > maxLenUnits)
          maxLenUnits = len;
        if ((len = strlen(qtyDesc[iQty].m_format.c_str())) > maxLenForms)
          maxLenForms = len;
        if ((len = strlen(qtyUCDs[iQty].c_str())) > maxLenUCDs)
          maxLenUCDs = len;
      }

      // Dump header
      Log(Log_2, "");
      Log(Log_2, "Catalogue descriptor:");
      Log(Log_2, "=====================");

      // Dump catalogue information
      Log(Log_2, " Catalogue input name .............: %s",
          in->inName.c_str());
      Log(Log_2, " Catalogue code ...................: %s",
          in->catCode.c_str());
      Log(Log_2, " Catalogue URL ....................: %s",
          in->catURL.c_str());
      Log(Log_2, " Catalogue name ...................: %s",
          in->catName.c_str());
      Log(Log_2, " Catalogue reference ..............: %s",
          in->catRef.c_str());
      Log(Log_2, " Catalogue table name .............: %s",
          in->tableName.c_str());
      Log(Log_2, " Catalogue table reference ........: %s",
          in->tableRef.c_str());
      Log(Log_2, " Number of objects in catalogue ...: %d",
          in->numTotal);
      Log(Log_2, " Number of loaded objects .........: %d",
          in->numLoad);
      Log(Log_2, " Number of quantities (columns) ...: %d",
          numQty);
      Log(Log_2, " Source ID column key .............: <%s>",
          in->col_id.c_str());
      switch (in->pos_type) {
      case Equatorial:
        Log(Log_2, " Position column keys (equatorial) : <%s> <%s>",
            in->col_ra.c_str(), in->col_dec.c_str());
        break;
      case Galactic:
        Log(Log_2, " Position column keys (galactic) ..: <%s> <%s>",
            in->col_glon.c_str(), in->col_glat.c_str());
        break;
      default:
        Log(Log_2, " Position column keys .............: "
                   "no position information found");
        break;
      }
      switch (in->col_e_type) {
      case Radius:
        Log(Log_2, " Position error column keys .......: <%s>",
            in->col_e_maj.c_str());
        break;
      case Ellipse:
        Log(Log_2, " Position error column keys .......: <%s> <%s> <%s>",
            in->col_e_maj.c_str(), in->col_e_min.c_str(), in->col_e_posang.c_str());
        break;
      case RaDec:
        Log(Log_2, " Position error column keys .......: <%s> <%s>",
            in->col_e_ra.c_str(), in->col_e_dec.c_str());
        break;
      case NoError:
      default:
        Log(Log_2, " Position error column keys .......: "
                   "no error information found");
        break;
      }
      switch (in->col_e_prob) {
      case Sigma_1:
        Log(Log_2, " Position error unit ..............: "
                   "1 sigma (68.269%%) (scale=%7.5f)", in->e_pos_scale);
        break;
      case Sigma_2:
        Log(Log_2, " Position error unit ..............: "
                   "2 sigma (95.450%%) (scale=%7.5f)", in->e_pos_scale);
        break;
      case Sigma_3:
        Log(Log_2, " Position error unit ..............: "
                   "3 sigma (99.730%%) (scale=%7.5f)", in->e_pos_scale);
        break;
      case Prob_68:
        Log(Log_2, " Position error unit ..............: 68%% (scale=%7.5f)",
            in->e_pos_scale);
        break;
      case Prob_90:
        Log(Log_2, " Position error unit ..............: 90%% (scale=%7.5f)",
            in->e_pos_scale);
        break;
      case Prob_95:
        Log(Log_2, " Position error unit ..............: 95%% (scale=%7.5f)",
            in->e_pos_scale);
        break;
      case Prob_99:
        Log(Log_2, " Position error unit ..............: 99%% (scale=%7.5f)",
            in->e_pos_scale);
        break;
      default:
        Log(Log_2, " Position error unit (default) ....: 95%% (scale=%7.5f)",
            in->e_pos_scale);
        break;
      }
      Log(Log_2, " 2D absolute 95%% position error ...: %7.5f deg", in->erposabs);

      // Dump information about catalogue quantitites
      if (par->logVerbose()) {
        for (iQty = 0; iQty < numQty; iQty++) {

          // Set quantity type
          switch (qtyTypes[iQty]) {
          case 0:
            qtyType = "vector";
            break;
          case 1:
            qtyType = "numerical";
            break;
          case 2:
            qtyType = "string";
            break;
          default:
            qtyType = "unknown";
            break;
          }

          // Dump quantity
          Log(Log_2,
              "  Quantity %3d ....................: %*s [%*s] (%*s) <%*s> (%s)",
              iQty+1,
              maxLenNames, qtyNames[iQty].c_str(),
              maxLenUnits, qtyUnits[iQty].c_str(),
              maxLenForms, qtyDesc[iQty].m_format.c_str(),
              maxLenUCDs,  qtyUCDs[iQty].c_str(),
              qtyType.c_str());

        } // endfor: looped over quantities
      } // endif: verbose level

    } while (0); // End of main do-loop

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Compute catalogue association posterior probability.
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] status Error status.
 * @param[in] quiet No information logging.
 *
 * To perform fast computation a sparse matrix is setup that hold the 
 * probabilities Pik'(H1|D). The rows of the matrix correspond to the NLAT
 * sources k. The columns of the matrix correspond to the Ncpt counterparts i.
 * The sparse matrix is stored in compressed sparse column format. The
 * following arrays are temporarily allocated to hold the sparse matrix
 * information:
 * tmp_prob   : holds all elements that are < 1, in column order
 * tmp_istart : holds the index of the first element of each column
 * tmp_k      : holds the row indices for each element.
 * Products over k (the LAT source indices) are quickly done by multiplication
 * over a column.
 ******************************************************************************/
Status Catalogue::compute_prob_post_cat(Parameters *par, Status status,
                                        int quiet) {

    // Temporary memory pointer
    double* tmp_prob   = NULL;  // Sparse matrix elements
    double* tmp_sum    = NULL;  // Normalization sum
    int*    tmp_k      = NULL;  // LAT source index k for each element
    int*    tmp_istart = NULL;  // First element index in each column

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::compute_prob_post_cat");

    // Single loop for common exit point
    do {

      // Dump header
      if (par->logNormal() && !quiet) {
        Log(Log_2, "");
        Log(Log_2, "Compute catalogue association probabilities:");
        Log(Log_2, "============================================");
      }

      // Allocate temporary memory for sparse matrix column start and
      // normalization sum
      tmp_istart = new int[m_cpt.numLoad+1];
      tmp_sum    = new double[m_cpt.numLoad];
      if (tmp_istart == NULL || tmp_sum == NULL) {
        status = STATUS_MEM_ALLOC;
        if (par->logTerse())
          Log(Error_2, "%d : Memory allocation failure.", (Status)status);
        continue;
      }

      // Initialise sparse matrix column start array
      for (int index = 0; index <= m_cpt.numLoad; ++index)
        tmp_istart[index] = 0;

      // Initialise normalization sums Si
      for (int index = 0; index < m_cpt.numLoad; ++index)
        tmp_sum[index] = 0.0;

      // Setup sparse matrix column start array. For this purpose we first
      // determine how often a specific index occurs. The result is stored
      // into the columns start array at position index+1. Then we build the
      // cumulative distribution by adding up successively the column start
      // values that appear before index. At the end the last elements
      // tmp_istart[m_cpt.numLoad] will contain the total number of counterparts
      // that occured.
      for (int k = 0; k < m_src.numLoad; ++k) {
        for (int i = 0; i < m_info[k].numRefine; ++i) {
          int index = m_info[k].cc[i].index + 1;
          tmp_istart[index] += 1;
        }
      }
      for (int index = 1; index <= m_cpt.numLoad; ++index)
        tmp_istart[index] += tmp_istart[index-1];

      // Determine the number of sparse matrix elements
      int num_elements = tmp_istart[m_cpt.numLoad];

      // Fall through if there are no counterparts
      if (num_elements < 1)
        continue;

      // Allocate memory for sparse matrix element values and indices
      tmp_prob   = new double[num_elements];
      tmp_k      = new int[num_elements];
      if (tmp_prob == NULL || tmp_k == NULL) {
        status = STATUS_MEM_ALLOC;
        if (par->logTerse())
          Log(Error_2, "%d : Memory allocation failure.", (Status)status);
        continue;
      }

      // Initialise element values (1.0) and indices (-1).
      // Indices of -1 signal that the value has not yet been set.
      for (int element = 0; element < num_elements; ++element) {
        tmp_prob[element] = 1.0;
        tmp_k[element]    = -1;
      }

      // Fill sparse matrix elements. For a given column we now search
      // for the first row that has an index of -1. This is then the next
      // free element. For security we check if all free elements are
      // exhausted. This should in principle never occur! If it occurs,
      // there is a logical error in this code.
      for (int k = 0; k < m_src.numLoad; ++k) {
        for (int i = 0; i < m_info[k].numRefine; ++i) {

          // Get sparse matrix column
          int col = m_info[k].cc[i].index;

          // Get start and stop element indices in sparse matrix
          int start = tmp_istart[col];
          int stop  = tmp_istart[col+1];

          // Set next free element
          int element = start;
          for ( ; element < stop; ++element) {
            if (tmp_k[element] == -1) {
              tmp_prob[element] = 1.0 - m_info[k].cc[i].prob_post_single;
              tmp_k[element]    = k;
              break;
            }
          }
          if (element >= stop) {
            status = STATUS_CAT_BAD_SPARSE;
            if (par->logTerse())
              Log(Error_2, "%d : No free element is sparse matrix.",
                  (Status)status);
            break;
          }

        }
        if (status != STATUS_OK)
          break;
      }
      if (status != STATUS_OK)
        continue;

      // Compute probability products. This is now done quickly by multiplying
      // down a column
      for (int k = 0; k < m_src.numLoad; ++k) {
        for (int i = 0; i < m_info[k].numRefine; ++i) {

          // Get sparse matrix column
          int col = m_info[k].cc[i].index;

          // Get start and stop element indices in sparse matrix
          int start = tmp_istart[col];
          int stop  = tmp_istart[col+1];

          // Compute products
          m_info[k].cc[i].prob_prod1 = 1.0; // all k'
          m_info[k].cc[i].prob_prod2 = 1.0; // all k' except of k
          for (int element = start; element < stop; ++element) {
            m_info[k].cc[i].prob_prod1   *= tmp_prob[element];
            if (k != tmp_k[element])
              m_info[k].cc[i].prob_prod2 *= tmp_prob[element];
          }

        }
      }

      // Initialise normalization sums with Pi(H-|D)
      for (int k = 0; k < m_src.numLoad; ++k) {
        for (int i = 0; i < m_info[k].numRefine; ++i) {
          int index      = m_info[k].cc[i].index;
          tmp_sum[index] = m_info[k].cc[i].prob_prod1;
        }
      }

      // Compute catalogue posterior probabilities and update normalization sum
      for (int k = 0; k < m_src.numLoad; ++k) {
        for (int i = 0; i < m_info[k].numRefine; ++i) {

          // Compute posterior probability
          m_info[k].cc[i].prob_post_cat = m_info[k].cc[i].prob_post_single *
                                          m_info[k].cc[i].prob_prod2;

          // Update normalization sum
          int index       = m_info[k].cc[i].index;
          tmp_sum[index] += m_info[k].cc[i].prob_post_cat;

        }
      }

      // Copy normalization sum into counterpart candidate field
      for (int k = 0; k < m_src.numLoad; ++k) {
        for (int i = 0; i < m_info[k].numRefine; ++i) {
          int index                 = m_info[k].cc[i].index;
          m_info[k].cc[i].prob_norm = tmp_sum[index];
        }
      }

      // Normalize catalogue posterior probabilities
      for (int k = 0; k < m_src.numLoad; ++k) {
        for (int i = 0; i < m_info[k].numRefine; ++i) {
          if (m_info[k].cc[i].prob_norm > 0.0)
            m_info[k].cc[i].prob_post_cat /= m_info[k].cc[i].prob_norm;
          else
            m_info[k].cc[i].prob_post_cat = 0.0;
        }
      }

      // Dump catalogue association probabilities
      if (par->logExplicit() && !quiet) {
        Log(Log_2,
            "  Source index   Counterpart index     P(H0|D) =>  P(Hk|D) "
            "    Pi_k' P  Pi_k'\\k P     Si");
        Log(Log_2,
            "  ------------  ------------------   --------- => ---------"
            "  ---------  ---------   -----");
        for (int k = 0; k < m_src.numLoad; ++k) {
          for (int i = 0; i < m_info[k].numRefine; ++i) {
            Log(Log_2,
                "  Source %5d: Counterpart %6d : %8.4f%% => %8.4f%%"
                " (%8.4f%%, %8.4f%%, %6.3f)",
                k+1, m_info[k].cc[i].index + 1,
                m_info[k].cc[i].prob_post_single*100.0,
                m_info[k].cc[i].prob_post_cat*100.0,
                m_info[k].cc[i].prob_prod1*100.0,
                m_info[k].cc[i].prob_prod2*100.0,
                m_info[k].cc[i].prob_norm);
          }
        }
      }

    } while (0); // End of main do-loop

    // Delete temporary memory
    if (tmp_prob   != NULL) delete [] tmp_prob;
    if (tmp_k      != NULL) delete [] tmp_k;
    if (tmp_istart != NULL) delete [] tmp_istart;
    if (tmp_sum    != NULL) delete [] tmp_sum;

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::compute_prob_post_cat (status=%d)",
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Compute unique catalogue association posterior probability.
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] status Error status.
 * @param[in] quiet No information logging.
 ******************************************************************************/
Status Catalogue::compute_prob_post(Parameters *par, Status status, int quiet) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::compute_prob_post");

    // Single loop for common exit point
    do {

      // Dump header
      if (par->logNormal() && !quiet) {
        Log(Log_2, "");
        Log(Log_2, "Compute unique catalogue association probabilities:");
        Log(Log_2, "===================================================");
      }

      // Initialise fraction of non unique sources
      m_fract_not_unique = 0.0;
      double num         = 0.0;

      // Initialise probability products for all counterpart candidates
      for (int k = 0; k < m_src.numLoad; ++k) {

        // Perform computations only if there are counterparts for this source
        if (m_info[k].numRefine > 0) {

          // Compute products
          for (int i = 0; i < m_info[k].numRefine; ++i) {
            m_info[k].cc[i].prob_prod1 = 1.0; // all i'
            m_info[k].cc[i].prob_prod2 = 1.0; // all i' except of i
            for (int ip = 0; ip < m_info[k].numRefine; ++ip) {
              if (i == ip)
                m_info[k].cc[i].prob_prod1  = (1.0 - m_info[k].cc[ip].prob_post_cat);
              else
                m_info[k].cc[i].prob_prod2 *= (1.0 - m_info[k].cc[ip].prob_post_cat);
            }
            m_info[k].cc[i].prob_prod1 *= m_info[k].cc[i].prob_prod2;
          }

          // Compute non-normalized posterior probabilities and normalization
          // factor
          double norm = m_info[k].cc[0].prob_prod1;
          for (int i = 0; i < m_info[k].numRefine; ++i) {
            m_info[k].cc[i].prob_post = m_info[k].cc[i].prob_post_cat *
                                        m_info[k].cc[i].prob_prod2;
            norm += m_info[k].cc[i].prob_post ;
          }

          // Compute unique association probabilities
          if (norm > 0.0) {
            for (int i = 0; i < m_info[k].numRefine; ++i) {
              m_info[k].cc[i].prob_post /= norm;
              m_info[k].cc[i].prob_norm  = norm;
            }
          }
          else {
            for (int i = 0; i < m_info[k].numRefine; ++i) {
              m_info[k].cc[i].prob_post = 0.0;
              m_info[k].cc[i].prob_norm = norm;
            }
          }

          // Update fraction of non unique sources
          m_fract_not_unique += (1.0 - norm);
          num                += 1.0;

        } // endif: there were counterparts

      } // endfor: looped over all sources

      // Compute fraction of non unique sources
      if (num > 0.0)
        m_fract_not_unique /= num;

      // Dump catalogue association probabilities
      if (par->logExplicit() && !quiet) {
        Log(Log_2,
            "  Source index   Counterpart index     P(Hk|D) =>  P(Hi|D) "
            "    Pi_i' P  Pi_i'\\i P     Sk");
        Log(Log_2,
            "  ------------  ------------------   --------- => ---------"
            "  ---------  ---------   -----");
        for (int k = 0; k < m_src.numLoad; ++k) {
          for (int i = 0; i < m_info[k].numRefine; ++i) {
            Log(Log_2,
                "  Source %5d: Counterpart %6d : %8.4f%% => %8.4f%%"
                " (%8.4f%%, %8.4f%%, %6.3f)",
                k+1, m_info[k].cc[i].index + 1,
                m_info[k].cc[i].prob_post_cat*100.0,
                m_info[k].cc[i].prob_post*100.0,
                m_info[k].cc[i].prob_prod1*100.0,
                m_info[k].cc[i].prob_prod2*100.0,
                m_info[k].cc[i].prob_norm);
          }
        }
      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::compute_prob_post (status=%d)",
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Compute association probability.
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] status Error status.
 ******************************************************************************/
Status Catalogue::compute_prob(Parameters *par, Status status) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::compute_prob");

    // Single loop for common exit point
    do {

      // Dump header
      if (par->logNormal()) {
        Log(Log_2, "");
        Log(Log_2, "Compute association probabilities:");
        Log(Log_2, "==================================");
      }

      // Initialise statistics
      m_sum_pid     = 0.0;
      m_sum_pc      = 0.0;
      m_sum_lr      = 0.0;
      m_sum_pid_thr = 0.0;
      m_sum_pc_thr  = 0.0;
      m_sum_lr_thr  = 0.0;
      m_num_claimed = 0.0;

      // Loop over all sources
      for (int k = 0; k < m_src.numLoad; ++k) {

        // Compute PROB
        status = cid_prob(par, &(m_info[k]), status);
        if (status != STATUS_OK) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to determine association probability.",
                (Status)status);
          break;
        }

        // Sort counterpart candidates by decreasing probability
        status = cid_sort(par, &(m_info[k]), m_info[k].numRefine, status);
        if (status != STATUS_OK) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to sort counterpart candidates.",
                (Status)status);
          break;
        }

        // Sum up probabilities before thresholding.
        for (int iCC = 0; iCC < m_info[k].numRefine; ++iCC) {
          m_sum_pid += m_info[k].cc[iCC].prob;
          m_sum_pc  += 1.0 - m_info[k].cc[iCC].prob;
          m_sum_lr  += m_info[k].cc[iCC].likrat;
        }

        // Determine the number of counterpart candidates above the probability
        // threshold
        int numUseCC = 0;
        for (int iCC = 0; iCC < m_info[k].numRefine; ++iCC) {
          if (m_info[k].cc[iCC].prob >= par->m_probThres)
            numUseCC++;
          else
            break;
        }

        // Apply the maximum number of counterpart threshold
        if (numUseCC > par->m_maxNumCpt)
          numUseCC = par->m_maxNumCpt;

        // Eliminate counterpart candidates below threshold.
        m_info[k].numClaimed = numUseCC;

        // Set unique counterpart candidate identifier (overwrite former ID)
        char cid[OUTCAT_MAX_STRING_LEN];
        for (int iCC = 0; iCC < m_info[k].numClaimed; ++iCC) {
          sprintf(cid, "CC_%5.5d_%5.5d", k+1, iCC+1);
          m_info[k].cc[iCC].id = cid;
        }

        // Sum up probabilities after thresholding
        for (int iCC = 0; iCC < m_info[k].numClaimed; ++iCC) {
          m_sum_pid_thr += m_info[k].cc[iCC].prob;
          m_sum_pc_thr  += 1.0 - m_info[k].cc[iCC].prob;
          m_sum_lr_thr  += m_info[k].cc[iCC].likrat;
          if (m_info[k].cc[iCC].likrat_div)
            m_num_lr_div += 1.0;
        }

        // Sum the total number of claimed associations
        m_num_claimed += double(m_info[k].numClaimed);

      } // endfor: looped over all sources

      // Compute reliability and completeness
      m_reliability  = (m_num_claimed > 0.0) ? m_sum_pid_thr/m_num_claimed : 0.0;
      m_completeness = (m_sum_pid     > 0.0) ? m_sum_pid_thr/m_sum_pid     : 0.0;

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::compute_prob (status=%d)",
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Perform catch-22 iterations
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] status Error status.
 ******************************************************************************/
Status Catalogue::catch22(Parameters *par, Status status) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::catch22");

    // Single loop for common exit point
    do {

      // Fall through if  catch-22 is not requested
      if (!par->m_catch22)
        continue;

      // Dump header
      if (par->logNormal()) {
        Log(Log_2, "");
        Log(Log_2, "Catch-22 computation of prior:");
        Log(Log_2, "==============================");
      }

      // Initialise convergence control
      double delta        = 0.0;
      double eps          = 0.0;
      double lambda       = 0.0;
      int    hit_boundary = 0;

      // Perform posterior probability iterations for catch-22 scheme
      for (m_iter = 0; m_iter < c_iter_max; ++m_iter) {

        // Break if prior run out of range
        if (hit_boundary)
          break;

        // Make next guess for prior
        double new_prior = 0.0;
        for (int k = 0; k < m_src.numLoad; ++k) {
          for (int i = 0; i < m_info[k].numRefine; ++i)
            new_prior += m_info[k].cc[i].prob_post;
        }
        new_prior /= double(m_cpt.numLoad);

        // Keep prior in range
        if (new_prior < m_prior_min) {
          new_prior    = m_prior_min;
          hit_boundary = 1;
        }
        if (new_prior > m_prior_max) {
          new_prior    = m_prior_max;
          hit_boundary = 1;
        }

        // Break if converged
        delta  = fabs(new_prior-m_prior);
        eps    = (m_prior > 0.0) ? delta/m_prior : 0.0;
        lambda = (m_prior > 0.0) ? new_prior/m_prior : 0.0;
        if (eps < 0.01 || delta < 1.0e-30) {
          if (par->logNormal()) {
            Log(Log_2, " Catch-22 converged prior prob. ...: %10.6f%%",
                m_prior*100.0);
            if (par->logExplicit()) {
              Log(Log_2, "  Prior change (Lambda) ...........: %10.3f", lambda);
              Log(Log_2, "  Relative convergence precision ..: %10.4f%%",
                  eps*100.0);
              Log(Log_2, "  Absolute convergence precision ..: %10.4f%%",
                  delta*100.0);
            }
          }
          break;
        }

        // Assign new prior
        m_prior = new_prior;

        // Dump new prior guess
        if (par->logNormal()) {
          if (par->logExplicit()) {
            Log(Log_2, " Iteration %4d ...................:", m_iter+1);
            Log(Log_2, "  New prior probability guess .....: %10.6f%%",
                m_prior*100.0);
            Log(Log_2, "  Prior change (Lambda) ...........: %10.3f", lambda);
          }
          else {
            Log(Log_2, " Iteration %4d ...................: %10.6f%%",
                m_iter+1, m_prior*100.0);
          }
        }

        // Re-compute PROB_POST_SINGLE for all sources
        for (int k = 0; k < m_src.numLoad; ++k) {

          // Update in-memory catalogue
          status = cfits_update(m_memFile, par, &(m_info[k]), m_info[k].numSelect,
                             status);
          if (status != STATUS_OK) {
            if (par->logTerse())
              Log(Error_2, "%d : Unable to update counterpart candidates"
                           " in-memory FITS catalogue.", (Status)status);
            break;
          }

          // Compute PROB_PRIOR
          status = cid_prob_prior(par, &(m_info[k]), status);
          if (status != STATUS_OK) {
            if (par->logTerse())
              Log(Error_2, "%d : Unable to re-compute PROB_PRIOR.",
                  (Status)status);
            break;
          }

          // Compute PROB_POST_SINGLE
          status = cid_prob_post_single(par, &(m_info[k]), status);
          if (status != STATUS_OK) {
            if (par->logTerse())
              Log(Error_2, "%d : Unable to re-compute PROB_POST_SINGLE.",
                  (Status)status);
            break;
          }

          // Copy over probabilities
          for (int i = 0; i < m_info[k].numSelect; ++i)
            m_info[k].cc[i].prob = m_info[k].cc[i].prob_post_single;

          // Sort probabilities
          status = cid_sort(par, &(m_info[k]), m_info[k].numSelect, status);
          if (status != STATUS_OK) {
            if (par->logTerse())
              Log(Error_2, "%d : Unable to sort probabilities.",
                  (Status)status);
            break;
          }

          // Select only relevant counterparts
          m_info[k].numRefine = 0;
          for (int i = 0; i < m_info[k].numSelect; ++i) {
            if (m_info[k].cc[i].prob_post_single <= c_prob_min)
              break;
            m_info[k].numRefine++;
          }

        } // endfor: looped over all sources
        if (status != STATUS_OK)
          break;

        // Compute probabilities for source catalogue association
        status = compute_prob_post_cat(par, status, 1);
        if (status != STATUS_OK) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to compute catalogue association"
                " probabilities.", (Status)status);
          break;
        }

        // Compute probabilities for unique source catalogue association
        status = compute_prob_post(par, status, 1);
        if (status != STATUS_OK) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to compute unique catalogue association"
                " probabilities.", (Status)status);
          break;
        }

      } // endfor: looped over posterior probability iterations

      // Signal if boundary was hit
      if (hit_boundary) {
        if (par->logNormal()) {
          Log(Warning_2, " Prior hit probability boundary ...: %10.6f%%"
              " [%10.6f%%, %10.6f%%]",
              m_prior*100.0, m_prior_min*100.0, m_prior_max*100.0);
        }
      }

      // Detect convergence problem
      if (m_iter >= c_iter_max) {
        if (par->logNormal()) {
          Log(Warning_2, " Catch-22 NON-CONVERGED prior .....: %10.6f%%",
              m_prior*100.0);
          if (par->logExplicit()) {
            Log(Warning_2, "  Prior change (Lambda) ...........: %10.3f", lambda);
            Log(Warning_2, "  Relative convergence precision ..: %10.4f%%",
                eps*100.0);
            Log(Warning_2, "  Absolute convergence precision ..: %10.4f%%",
                delta*100.0);
          }
        }
      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::catch22 (status=%d)",
          status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Dump counterpart identification results
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] status Error status.
 *
 * The results are stored in
 * the array m_cpt_stat[m_src.numLoad * (m_num_sel+1)]
 * the vector m_src_cpts[m_src.numLoad] and
 * the vector m_cpt_names[m_src.numLoad].
 * m_cpt_stat contains for each source a vector that traces the number of
 * counterparts that survive a given step. The first elements is the number
 * of counterpart candidates that came out of the filter step. The following
 * elements are the number of counterpart that survived the various selection
 * criteria. m_src_cpts contains for each source the number of counterparts
 * that survived the refine step. The names of these counterparts and their
 * associated probabilities are stored in the vector m_cpt_names.
 ******************************************************************************/
Status Catalogue::dump_results(Parameters *par, Status status) {

    // Declare local variables
    ObjectInfo *src;

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::dump_results");

    // Single loop for common exit point
    do {

      // Fall through in case of an error
      if (status != STATUS_OK)
        continue;

      // Dump header
      Log(Log_2, "");
      Log(Log_2, "Counterpart association results:");
      Log(Log_2, "================================");

      // Build header
      char add[256];
      char select[256] = "";
      for (int iSel = 0; iSel < m_num_Sel; ++iSel) {
        sprintf(add, " Sel%2.2d", iSel+1);
        strcat(select, add);
      }
      sprintf(add, " Select");
      strcat(select, add);
      sprintf(add, " Refine");
      strcat(select, add);
      sprintf(add, "  Claim");
      strcat(select, add);
      sprintf(add, "  Final");
      strcat(select, add);

      // Dump header
      Log(Log_2, "                                      Filter%s", select);

      // Loop over all sources
      int n_src_assoc = 0;
      for (int iSrc = 0; iSrc < m_src.numLoad; ++iSrc) {

        // Get pointer to source object
        src = &(m_src.object[iSrc]);

        // Build selection string
        sprintf(select, " %6d", m_info[iSrc].numFilter);
        for (int iSel = 0; iSel < m_num_Sel; ++iSel) {
          sprintf(add, " %5d", m_cpt_stat[iSrc*(m_num_Sel+1) + iSel+1]);
          strcat(select, add);
        }
        sprintf(add, " %6d", m_info[iSrc].numSelect);
        strcat(select, add);
        sprintf(add, " %6d", m_info[iSrc].numRefine);
        strcat(select, add);
        sprintf(add, " %6d", m_info[iSrc].numClaimed);
        strcat(select, add);
        sprintf(add, " %6d", m_info[iSrc].numFinalSel);
        strcat(select, add);

        // Dump information
        Log(Log_2, " Source %5d %20s : %s %s",
            iSrc+1, src->name.c_str(), select, m_cpt_names[iSrc].c_str());

        // Collect number of associated sources
        if (m_info[iSrc].numFinalSel > 0)
          n_src_assoc++;

      } // endfor: looped over all sources

      // Compute some results
      double f_assoc  = (m_src.numLoad > 0) ? double(n_src_assoc) /
                                              double(m_src.numLoad) : 0.0;

      // Dump summary
      Log(Log_2, "");
      Log(Log_2, "Counterpart association summary:");
      Log(Log_2, "================================");
      Log(Log_2, " Number of sources .............................: %10d",
          m_src.numLoad);
      Log(Log_2, " Number of sources with at least one counterpart: %10d (%.3f%%)",
          n_src_assoc, f_assoc*100.0);
      Log(Log_2, " Number of claimed identifications .............: %10d",
          int(m_num_claimed));
      Log(Log_2, " Expected number of true associations ..........: %10.1f (sum of posterior probabilities)",
          m_sum_pid);
      Log(Log_2, " Expected number of correct identifications ....: %10.1f (sum of posterior probabilities above threshold)",
          m_sum_pid_thr);
      Log(Log_2, " Expected number of spurious identifications ...: %10.1f (sum of 1-posterior probabilitis above threshold)",
          m_sum_pc_thr);
      Log(Log_2, " Expected number of missed identifications .....: %10.1f (difference between total posterior sum and sum above threshold)",
          m_sum_pid-m_sum_pid_thr);
      Log(Log_2, " Expected number of non-unique counterparts ....: %10.1f (%.3f%%)",
          m_fract_not_unique*m_num_claimed, m_fract_not_unique*100.0);
      Log(Log_2, " Reliability of identifications ................: %10.3f%%",
          m_reliability*100.0);
      Log(Log_2, " Completeness of identifications ...............: %10.3f%%",
          m_completeness*100.0);
      Log(Log_2, " Total log likelihood-ratio ....................: %10.3f"
                 " (before threshold: %.3f)", m_sum_lr_thr, m_sum_lr);
      Log(Log_2, " Number of associations with divergent LR ......: %10d",
          int(m_num_lr_div));
      if (par->m_catch22) {
        Log(Log_2, " Catch-22 converged prior probability ..........: %10.6f%%",
            m_prior*100.0);
        Log(Log_2, " Number of Catch-22 iterations .................: %10.d",
            m_iter);
      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::dump_results (status=%d)", status);

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Build counterpart catalogue
 *
 * @param[in] par Pointer to gtsrcid parameters.
 * @param[in] status Error status.
 *
 * Main driver method that performs counterpart associations.
 *
 * \verbatim
 * build
 *   |
 *   +-- get_input_descriptor (get source catalogue input descriptior)
 *   |
 *   +-- get_input_descriptor (get counterpart catalogue input descriptior)
 *   |
 *   +-- cfits_create (create output catalogue)
 *   |
 *   +-- get_input_catalogue (get source catalogue)
 *   |
 *   N-- cid_single (perform single source association)
 *   |   |
 *   |   +-- cid_filter (filter step)
 *   |   |   |
 *   |   |   +-- get_input_catalogue (get counterpart catalogue)
 *   |   |
 *   |   +-- cid_select (select counterparts)
 *   |   |   |
 *   |   |   +-- cfits_select (select output catalogue entries)
 *   |   |
 *   |   +-- cid_refine (refine step)
 *   |       |
 *   |       +-- cfits_clear (clear in-memory catalogue)
 *   |       |
 *   |       +-- cfits_add (add quantities to in-memory catalogue)
 *   |       |
 *   |       +-- cid_prob_pos (compute PROB_POS & PDF_POS)
 *   |       |
 *   |       +-- cid_prob_prior (compute PROB_PRIOR)
 *   |       |
 *   |       +-- cid_prob_chance (compute PROB_CHANCE & PDF_CHANCE)
 *   |       |
 *   |       +-- cid_prob_post_single (compute PROB_POST)
 *   |
 *   +-- compute_prob_post_cat (compute catalogue association probabilities)
 *   |
 *   +-- compute_prob_post (compute unique catalogue association probabilities)
 *   |
 *   +-- catch22 (perform Catch22 iterations)
 *   |
 *   +-- compute_prob (compute association probability)
 *   |
 *   N-- cfits_add (add counterpart candidates to output catalogue)
 *   |
 *   +-- cfits_eval (evaluate output catalogue quantities)
 *   |
 *   +-- cfits_select (select output catalogue entries)
 *   |
 *   +-- cfits_collect (collect counterpart results)
 *   |
 *   +-- cfits_set_pars (set run parameter keywords)
 *   |
 *   +-- cfits_save (save output catalogue)
 *   |
 *   +-- dump_results (dump results)
 * \endverbatim
 ******************************************************************************/
Status Catalogue::build(Parameters *par, Status status) {

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " ==> ENTRY: Catalogue::build");

    // Single loop for common exit point
    do {

      // Fall through in case of an error
      if (status != STATUS_OK)
        continue;

      // Dump header
      if (par->logNormal()) {
        Log(Log_2, "");
        Log(Log_2, "Prepare catalogues:");
        Log(Log_2, "===================");
      }

      // Get input catalogue descriptors
      status = get_input_descriptor(par, par->m_srcCatName, &m_src, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to load source catalogue '%s' descriptor.",
              (Status)status, par->m_srcCatName.c_str());
        continue;
      }
      status = get_input_descriptor(par, par->m_cptCatName, &m_cpt, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to load counterpart catalogue '%s'"
              " descriptor.", (Status)status, par->m_cptCatName.c_str());
        continue;
      }
      if (par->logTerse()) {
        status = dump_descriptor(par, &m_src, status);
        if (status != STATUS_OK)
          continue;
      }
      if (par->logTerse()) {
        status = dump_descriptor(par, &m_cpt, status);
        if (status != STATUS_OK)
          continue;
      }

      // Optionally read counterpart density map
      if (par->m_cptDensFile.length() > 0) {
        m_has_density = 0;
        if (par->logNormal()) {
          Log(Log_2, "");
          Log(Log_2, "Read counterpart catalogue density map:");
          Log(Log_2, "=======================================");
        }
        try {
          m_density     = GHealpix(par->m_cptDensFile);
          m_has_density = 1;
          if (par->logNormal()) {
            Log(Log_2, " Filename .........................: %s", par->m_cptDensFile.c_str());
            Log(Log_2, " Nside (number of divisions) ......: %d", m_density.nside());
          }
        }
        catch (std::exception &e) {
          if (par->logTerse()) {
            Log(Warning_3, "Error occured while loading file '%s' (standard exception)",
                            par->m_cptDensFile.c_str());
          }
        }
        catch (std::string str) {
          if (par->logTerse()) {
            Log(Warning_3, "Error occured while loading file '%s': %s.",
                            par->m_cptDensFile.c_str(), str.c_str());
          }
        }
        catch (...) {
          if (par->logTerse()) {
            Log(Warning_3, "Error occured while loading file '%s' (unknown exception)",
                            par->m_cptDensFile.c_str());
          }
        }
      }

      // Create FITS catalogue in memory
      status = cfits_create(&m_memFile, "mem://gtsrcid", par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to create FITS memory catalogue"
              " 'mem://gtsrcid'.", (Status)status);
        continue;
      }

      // Create FITS output catalogue on disk
      status = cfits_create(&m_outFile, (char*)par->m_outCatName.c_str(), par,
                            status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to create FITS output catalogue '%s'.",
              (Status)status, par->m_outCatName.c_str());
        continue;
      }

      // Dump header
      if (par->logNormal()) {
        Log(Log_2, "");
        Log(Log_2, "Build counterpart candidate catalogue:");
        Log(Log_2, "======================================");
      }

      // Load source catalogue
      status = get_input_catalogue(par, &m_src, par->m_srcPosError, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to load source catalogue '%s' data.",
              (Status)status, par->m_srcCatName.c_str());
        continue;
      }
      else {
        if (par->logVerbose())
          Log(Log_2, " Source catalogue loaded.");
      }

      // Stop if the source catalogue is empty
      if (m_src.numLoad < 1) {
        status = STATUS_CAT_EMPTY;
        if (par->logTerse())
          Log(Error_2, "%d : Source catalogue is empty. Stop", (Status)status);
        continue;
      }
      else {
        if (par->logVerbose())
          Log(Log_2, " Source catalogue contains %d sources.", m_src.numLoad);
      }

      // Load counterpart catalogue
      status = get_input_catalogue(par, &m_cpt, par->m_cptPosError, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to load counterpart catalogue '%s' data.",
              (Status)status, par->m_cptCatName.c_str());
        continue;
      }
      else {
        if (par->logVerbose())
          Log(Log_2, " Counterpart catalogue loaded.");
      }

      // Stop if the counterpart catalogue is empty
      if (m_cpt.numLoad < 1) {
        status = STATUS_CAT_EMPTY;
        if (par->logTerse())
          Log(Error_2, "%d : Counterpart catalogue is empty. Stop.", (Status)status);
        continue;
      }
      else {
        if (par->logVerbose())
          Log(Log_2, " Counterpart catalogue contains %d sources.", m_cpt.numLoad);
      }

      // Allocate source information
      m_info = new SourceInfo[m_src.numLoad];
      if (m_info == NULL) {
        status = STATUS_MEM_ALLOC;
        if (par->logTerse())
          Log(Error_2, "%d : Memory allocation failure.", (Status)status);
        continue;
      }

      // Initialise source information
      for (int iSrc = 0; iSrc < m_src.numLoad; ++iSrc) {
        m_info[iSrc].iSrc         = iSrc;
        m_info[iSrc].info         = &(m_src.object[iSrc]);
        m_info[iSrc].numFilter    = 0;
        m_info[iSrc].numSelect    = 0;
        m_info[iSrc].numRefine    = 0;
        m_info[iSrc].numClaimed   = 0;
        m_info[iSrc].numFinalSel  = 0;
        m_info[iSrc].cc           = NULL;
        m_info[iSrc].filter_rad   = 0.0;
        m_info[iSrc].ring_rad_min = 0.0;
        m_info[iSrc].ring_rad_max = 0.0;
        m_info[iSrc].omega        = 0.0;
      }

      // Allocate list of selected counterparts
      m_cpt_sel = new int[m_cpt.numLoad];
      if (m_cpt_sel == NULL) {
        status = STATUS_MEM_ALLOC;
        if (par->logTerse())
          Log(Error_2, "%d : Memory allocation failure.", (Status)status);
        continue;
      }

      // Determine number of quantity selection criteria
      m_num_Sel = par->m_select.size();

      // Set vectors dimensions
      m_cpt_names = std::vector<std::string>(m_src.numLoad);

      // Allocate selection statistics
      m_cpt_stat = new int[m_src.numLoad*(m_num_Sel+1)];
      if (m_cpt_stat == NULL) {
        status = STATUS_MEM_ALLOC;
        if (par->logTerse())
          Log(Error_2, "%d : Memory allocation failure.", (Status)status);
        continue;
      }
      for (int iSrc = 0; iSrc < m_src.numLoad; ++iSrc) {
        for (int iSel = 0; iSel <= m_num_Sel; ++iSel)
          m_cpt_stat[iSrc*(m_num_Sel+1) + iSel] = 0;
      }

      // Get plausible counterpart candidates and compute PROB_POST_SINGLE
      // for them
      for (int iSrc = 0; iSrc < m_src.numLoad; ++iSrc)
        status = cid_source(par, &(m_info[iSrc]), status);
      if (status != STATUS_OK)
        continue;

      // Compute probabilities for source catalogue association
      status = compute_prob_post_cat(par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to compute catalogue association"
              " probabilities.", (Status)status);
        continue;
      }

      // Compute probabilities for unique source catalogue association
      status = compute_prob_post(par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to compute unique catalogue association"
              " probabilities.", (Status)status);
        continue;
      }

      // Perform catch-22 iterations
      status = catch22(par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to perform catch-22 iterations.",
              (Status)status);
        continue;
      }

      // Compute association probabilities
      status = compute_prob(par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to compute association probabilities.",
              (Status)status);
        continue;
      }

      // Save claimed counterpart candidates for all sources
      for (int iSrc = 0; iSrc < m_src.numLoad; ++iSrc) {
        status = cfits_add(m_outFile, par, &(m_info[iSrc]), 
                           m_info[iSrc].numClaimed, status);
        if (status != STATUS_OK) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to add counterpart candidates for source"
                " %d to FITS output catalogue '%s'.",
                (Status)status, iSrc+1, par->m_outCatName.c_str());
          break;
        }
      }
      if (status != STATUS_OK)
        continue;

      // Collect statistics (used to build counterpart names)
//      std::vector<int> stat;
//      status = cfits_collect(m_outFile, par, stat, status);
//      if (status != STATUS_OK) {
//        if (par->logTerse())
//          Log(Error_2, "%d : Unable to collect statistics from catalogue.",
//                        (Status)status);
//        continue;
//      }

      // Evaluate output catalogue quantities
      if (par->logNormal()) {
        Log(Log_2, "");
        Log(Log_2, "Evaluate new output catalogue quantities:");
        Log(Log_2, "=========================================");
      }
      status = cfits_eval(m_outFile, par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to evaluate output catalogue quantities.",
              (Status)status);
        continue;
      }

      // Perform final selection
      if (par->logNormal()) {
        Log(Log_2, "");
        Log(Log_2, "Select output catalogue sources:");
        Log(Log_2, "================================");
      }
      status = cfits_select(m_outFile, par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to select output catalogue sources.",
              (Status)status);
        continue;
      }

      // Collect statistics (used to build counterpart names)
      std::vector<int> stat;
      status = cfits_collect(m_outFile, par, stat, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to collect statistics from catalogue.",
                        (Status)status);
        continue;
      }

      // Determine number of finally selected sources
      for (int iSrc = 0; iSrc < m_src.numLoad; ++iSrc)
        m_info[iSrc].numFinalSel  = stat[iSrc];

      // Write parameters as keywords to catalogue
      status = cfits_set_pars(m_outFile, par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to write parameters to output catalogue.",
              (Status)status);
        continue;
      }

      // Save output catalogue counterparts
      status = cfits_save(m_outFile, par, status);
      if (status != STATUS_OK) {
        if (par->logTerse())
          Log(Error_2, "%d : Unable to save output catalogue.",
              (Status)status);
        continue;
      }

      // Dump counterpart results
      if (par->logTerse()) {
        status = dump_results(par, status);
        if (status != STATUS_OK) {
          if (par->logTerse())
            Log(Error_2, "%d : Unable to dump counterpart results.",
                (Status)status);
          continue;
        }
      }

    } while (0); // End of main do-loop

    // Debug mode: Entry
    if (par->logDebug())
      Log(Log_0, " <== EXIT: Catalogue::build (status=%d)", status);

    // Return status
    return status;

}


/* Namespace ends ___________________________________________________________ */
}
