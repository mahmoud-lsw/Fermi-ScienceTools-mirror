/*------------------------------------------------------------------------------
Id ........: $Id: Parameters.cxx,v 1.17 2008/07/08 18:43:15 jurgen Exp $
Author ....: $Author: jurgen $
Revision ..: $Revision: 1.17 $
Date ......: $Date: 2008/07/08 18:43:15 $
--------------------------------------------------------------------------------
$Log: Parameters.cxx,v $
Revision 1.17  2008/07/08 18:43:15  jurgen
Remove GtApp, parametrize prefix symbol and update unit test1

Revision 1.16  2008/04/24 14:55:17  jurgen
Implement simple FoM scheme

Revision 1.15  2008/04/23 15:09:21  jurgen
Check for catch-22 in upper case

Revision 1.14  2008/04/18 20:50:33  jurgen
Implement catch-22 scheme for prior probability calculation and compute log likelihood-ratio instead of likelihood ratio (avoid numerical problems)

Revision 1.13  2008/03/26 13:37:10  jurgen
Generalize probability calculation and implement Bayesian method

Revision 1.12  2008/03/21 09:10:12  jurgen
Enhance code documentation.

Revision 1.11  2008/03/20 21:56:26  jurgen
implement local counterpart density

Revision 1.10  2007/10/09 08:17:40  jurgen
Correctly interpret positional errors and correctly evaluate PROB_POS
as likelihood

Revision 1.9  2007/10/08 11:02:25  jurgen
Implement search for catalogue table information and handle different
position error types

Revision 1.8  2007/10/02 22:01:16  jurgen
Change parameter name maxNumCtp to maxNumCpt

Revision 1.7  2006/02/09 15:51:42  jurgen
Remove unreferenced local variable 'pos'

Revision 1.6  2006/02/07 11:10:51  jurgen
Suppress catalogAccess verbosity

Revision 1.5  2006/02/06 13:26:14  jurgen
Remove whitespace in input parameter strings

Revision 1.4  2006/02/03 12:14:52  jurgen
New version that allows additional probabilities to be taken
into account. The code has been considerably reorganised. Also
catalogue column prefixes are now handled differently.

Revision 1.3  2006/02/01 13:33:36  jurgen
Tried to fix Win32 compilation bugs.
Change revision number to 1.3.2.
Replace header information with CVS typeset information.

------------------------------------------------------------------------------*/
/**
 * @file Parameters.h
 * @brief Parameters class implementation.
 * @author J. Knodlseder
 */

/* Includes _________________________________________________________________ */
#include <stdio.h>                       // for "sprintf" function
#include "sourceIdentify.h"
#include "Parameters.h"
#include "Log.h"                         // for parameter dumping/errors
#include "Catalogue.h"

/* Namespace definition _____________________________________________________ */
namespace sourceIdentify {


/* Globals __________________________________________________________________ */
int g_u9_verbosity;


/* Type defintions __________________________________________________________ */


/* Prototypes _______________________________________________________________ */
std::string trim(std::string str);


/*============================================================================*/
/*                              Private functions                             */
/*============================================================================*/

/**************************************************************************//**
 * @brief Remove leading and trailing whitespace from string
 *
 * @param[in] str String from which whitespace is to be removed.
 ******************************************************************************/
std::string trim(std::string str) {

    // Declare variables
    char const* delims = " \t\r\n";

    // Declare variables
    std::string::size_type notwhite;

    // Trim leading whitespace
    notwhite = str.find_first_not_of(delims);
    str.erase(0,notwhite);

    // Trim trailing whitespace
    notwhite = str.find_last_not_of(delims);
    str.erase(notwhite+1);

    // Return string
    return str;

}


/*============================================================================*/
/*                          Low-level parameter methods                       */
/*============================================================================*/

/**************************************************************************//**
 * @brief Initialise class memory
 ******************************************************************************/
void Parameters::init_memory(void) {

    // Declare local variables

    // Single loop for common exit point
    do {

      // Intialise private members
      m_srcCatName.clear();
      m_srcCatPrefix.clear();
      m_srcCatQty.clear();
      m_cptCatName.clear();
      m_cptCatPrefix.clear();
      m_cptCatQty.clear();
      m_cptDensFile.clear();
      m_outCatName.clear();
      m_outCatQtyName.clear();
      m_outCatQtyFormula.clear();
      m_probMethod.clear();
      m_probPrior.clear();
      m_select.clear();
      m_probThres   = 0.0;
      m_srcPosError = 0.0;
      m_cptPosError = 0.0;
      m_maxNumCpt   = 0;
      m_catch22     = 0;
      m_chatter     = 0;
      m_clobber     = 0;
      m_debug       = 0;
      m_mode.clear();

    } while (0); // End of main do-loop

    // Return
    return;

}


/**************************************************************************//**
 * @brief Free class memory
 ******************************************************************************/
void Parameters::free_memory(void) {

    // Declare local variables

    // Single loop for common exit point
    do {

      // Initialise memory
      init_memory();

    } while (0); // End of main do-loop

    // Return
    return;

}


/**************************************************************************//**
 * @brief Load parameters from task parameter file
 *
 * @param[in] pars Task parameters.
 * @param[in] status Error status.
 *
 * The probability method string is split into its elements (separated by '*').
 * The following probability names are detected:
 * 'POSITION':  Probability based on position, using Gaussian errors.
 * 'POS-GAUSS': Same as 'POSITION'.
 * 'POS-EXP':   Probability based on position, using exponential errors.
 * 'CHANCE':    Inclusion of local chance coincidence probability.
 * All other names correspond to columns in the catalogues.
 ******************************************************************************/
Status Parameters::load(st_app::AppParGroup &pars, Status status) {

    // Declare local variables
    char                   parname[MAX_CHAR];
    std::string::size_type len;
    std::string::size_type pos;
    std::string::size_type len_name;
    std::string::size_type start_formula;
    std::string::size_type len_formula;

    // Single loop for common exit point
    do {

      // Fall through in case of an error
      if (status != STATUS_OK)
        continue;

      // Prompt for parameters
      pars.Prompt();

      // Save parameters
      pars.Save();

      // Recover task parameters
      std::string s_srcCatName   = pars["srcCatName"];
      std::string s_srcCatPrefix = pars["srcCatPrefix"];
      std::string s_srcCatQty    = pars["srcCatQty"];
      std::string s_cptCatName   = pars["cptCatName"];
      std::string s_cptCatPrefix = pars["cptCatPrefix"];
      std::string s_cptCatQty    = pars["cptCatQty"];
      std::string s_cptDensFile  = pars["cptDensFile"];
      std::string s_outCatName   = pars["outCatName"];
      std::string s_probMethod   = pars["probMethod"];
      std::string s_probPrior    = pars["probPrior"];
      std::string s_FoM          = pars["fom"];
      std::string s_mode         = pars["mode"];
      m_srcCatName               = trim(s_srcCatName);
      m_srcCatPrefix             = OUTCAT_PRE_STRING + s_srcCatPrefix + "_";
      m_srcCatQty                = s_srcCatQty;
      m_srcPosError              = pars["srcPosError"];
      m_cptCatName               = trim(s_cptCatName);
      m_cptCatPrefix             = OUTCAT_PRE_STRING + s_cptCatPrefix + "_";
      m_cptCatQty                = s_cptCatQty;
      m_cptPosError              = pars["cptPosError"];
      m_cptDensFile              = trim(s_cptDensFile);
      m_outCatName               = trim(s_outCatName);
      m_probMethod               = trim(s_probMethod);
      m_probPrior                = trim(s_probPrior);
      m_FoM                      = trim(s_FoM);
      m_probThres                = pars["probThres"];
      m_maxNumCpt                = pars["maxNumCpt"];
      m_chatter                  = pars["chatter"];
      m_clobber                  = pars["clobber"];
      m_debug                    = pars["debug"];
      m_mode                     = s_mode;

      // Set U9 verbosity
      if (m_debug)
        g_u9_verbosity = 3;
      else
        g_u9_verbosity = 0;

      // Retrieve new output quantities and decompose them into quantity name
      // and evaluation string
      for (int i = MIN_OUTCAT_QTY; i <= MAX_OUTCAT_QTY; ++i) {

        // Extract parameter name
        sprintf(parname, "outCatQty%2.2d", i);
        std::string outCatQty = pars[parname];

        // Fall through if parameter is empty
        outCatQty = trim(outCatQty);
        len       = outCatQty.length();
        if (len < 1)
          continue;

        // Decompose string in part before and after "=" symbol
        pos           = outCatQty.find("=",0);
        len_name      = pos;
        start_formula = pos + 1;
        len_formula   = len - start_formula;

        // Catch invalid parameters
        if (pos == std::string::npos) {
          status = STATUS_PAR_BAD_PARAMETER;
          Log(Error_2, "%d : No equality symbol found in new output catalogue"
              " quantity string <%s='%s'>.", 
              (Status)status, parname, outCatQty.c_str());
          break;
        }
        if (len_name < 1) {
          status = STATUS_PAR_BAD_PARAMETER;
          Log(Error_2, "%d : No quantity name found for new output catalogue"
              " quantity <%s='%s'>.", 
              (Status)status, parname, outCatQty.c_str());
          break;
        }
        if (len_formula < 1) {
          status = STATUS_PAR_BAD_PARAMETER;
          Log(Error_2, "%d : No quantity evaluation string found for new"
              " output catalogue quantity <%s='%s'>.", 
              (Status)status, parname, outCatQty.c_str());
          break;
        }

        // Set name and formula (remove whitespace)
        m_outCatQtyName.push_back(trim(outCatQty.substr(0, len_name)));
        m_outCatQtyFormula.push_back(trim((outCatQty.substr(start_formula,
                                                            len_formula))));

      } // endfor: looped over quantities
      if (status != STATUS_OK)
        continue;

      // Retrieve selection strings
      for (int i = MIN_OUTCAT_SEL; i <= MAX_OUTCAT_SEL; ++i) {
        sprintf(parname, "select%2.2d", i);
        std::string select = pars[parname];
        select             = trim(select);
        len                = select.length();
        if (len > 0) {
          m_select.push_back(select);
        }
      }
      if (status != STATUS_OK)
        continue;

      // Check for catch-22
      std::string u_probPrior = upper(s_probPrior);
      if ((u_probPrior.find("CATCH-22",0) != std::string::npos) ||
          (u_probPrior.find("CATCH22",0)  != std::string::npos))
        m_catch22 = 1;

    } while (0); // End of main do-loop

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Dump task parameters into log file
 *
 * @param[in] status Error status.
 ******************************************************************************/
Status Parameters::dump(Status status) {

    // Declare local variables
    std::string::size_type i;
    std::string::size_type n;

    // Single loop for common exit point
    do {

      // Dump task parameters
      Log(Log_1, "Task Parameters:");
      Log(Log_1, "================");
      Log(Log_1, " Source catalogue filename ........: %s", m_srcCatName.c_str());
      Log(Log_1, " Source catalogue prefix ..........: %s", m_srcCatPrefix.c_str());
      Log(Log_1, " Source catalogue quantities ......: %s", m_srcCatQty.c_str());
      Log(Log_1, " Source catalogue uncertainty .....: %.3f arcmin",
          m_srcPosError*60.0);
      Log(Log_1, " Counterpart catalogue name .......: %s", m_cptCatName.c_str());
      Log(Log_1, " Counterpart catalogue prefix .....: %s", m_cptCatPrefix.c_str());
      Log(Log_1, " Counterpart catalogue quantities .: %s", m_cptCatQty.c_str());
      Log(Log_1, " Counterpart catalogue uncertainty : %.3f arcmin",
          m_cptPosError*60.0);
      if (m_cptDensFile.length() > 0)
        Log(Log_1, " Counterpart catalogue density file: %s", m_cptDensFile.c_str());
      else
        Log(Warning_1, " Counterpart catalogue density file: not used");
      Log(Log_1, " Output catalogue name ............: %s", m_outCatName.c_str());
      Log(Log_1, " Association probability ..........: PROB = %s",
          m_probMethod.c_str());
      if (m_catch22)
        Log(Log_1, " Counterpart association prior ....: CATCH-22");
      else
        Log(Log_1, " Counterpart association prior ....: PROB_PRIOR = %s",
            m_probPrior.c_str());
      Log(Log_1, " Probability threshold ............: %.3e", m_probThres);
      Log(Log_1, " Max. number of cpts per source ...: %d", m_maxNumCpt);
      if (m_FoM.length() > 0)
        Log(Log_1, " Figure of merit ..................: FoM = %s", m_FoM.c_str());
      else
        Log(Warning_1, " Figure of merit ..................: not used");
      if ((n = m_outCatQtyName.size()) > 0) {
        for (i = 0; i < n; ++i) {
          Log(Log_1, " New output catalogue quantity %2d .: %s = %s",
              i+1, m_outCatQtyName[i].c_str(), m_outCatQtyFormula[i].c_str());
        }
      }
      if ((n = m_select.size()) > 0) {
        for (i = 0; i < n; i++) {
          Log(Log_1, " Output catalogue selection %2d ....: %s",
              i+1, m_select[i].c_str());
        }
      }
      Log(Log_1, " Chatter level of output ..........: %d", m_chatter);
      Log(Log_1, " U9 verbosity .....................: %d", g_u9_verbosity);
      Log(Log_1, " Clobber ..........................: %d", m_clobber);
      Log(Log_1, " Debugging mode activated .........: %d", m_debug);
      Log(Log_1, " Mode of automatic parameters .....: %s", m_mode.c_str());

    } while (0); // End of main do-loop

    // Return status
    return status;

}


/* Namespace ends ___________________________________________________________ */
}
