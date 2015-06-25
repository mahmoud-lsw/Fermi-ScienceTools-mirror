/*------------------------------------------------------------------------------
Id ........: $Id: Log.cxx,v 1.6 2012/01/17 11:07:42 jurgen Exp $
Author ....: $Author: jurgen $
Revision ..: $Revision: 1.6 $
Date ......: $Date: 2012/01/17 11:07:42 $
--------------------------------------------------------------------------------
$Log: Log.cxx,v $
Revision 1.6  2012/01/17 11:07:42  jurgen
Correct format string for logging
Change CESR to IRAP
Update version number to v2r3p4

Revision 1.5  2008/03/21 09:10:12  jurgen
Enhance code documentation.

Revision 1.4  2007/09/21 14:29:03  jurgen
Correct memory bug and updated test script

Revision 1.3  2006/02/01 13:33:36  jurgen
Tried to fix Win32 compilation bugs.
Change revision number to 1.3.2.
Replace header information with CVS typeset information.

------------------------------------------------------------------------------*/
/**
 * @file Log.cxx
 * @brief Logging interface implementation.
 * @author J. Knodlseder
 */

/* Includes _________________________________________________________________ */
#include <stdio.h>      // for "FILE" type
#include <string.h>     // for "memcpy" function
#include <stdarg.h>     // for "va_list" type
#include <time.h>       // for time functions
#include "sourceIdentify.h"
#include "Log.h"


/* Namespace definition _____________________________________________________ */
namespace sourceIdentify {


/* Globals __________________________________________________________________ */
char  gLogTaskName[100];
FILE *gLogFilePtr = NULL;


/* Type defintions __________________________________________________________ */


/* Prototypes _______________________________________________________________ */


/**************************************************************************//**
 * @brief Initialise task logging
 *
 * @param[in] logName Name of log file.
 * @param[in] taskName Name of task.
 * @param[in] status Error status.
 ******************************************************************************/
Status LogInit(const char *logName, const char *taskName, Status status) {

    // Declare (and initialise) variables

    // Main do-loop to fall through in case of an error
    do {

      // Clean log task name
      sprintf(gLogTaskName, "%s", "");

      // If there is already a log file opened then close it first
      if (gLogFilePtr != NULL) {
        status = LogClose(status);
        if (status != STATUS_OK)
          continue;
        else
          gLogFilePtr = NULL;
      }

      // Open log file
      if (gLogFilePtr == NULL)
        gLogFilePtr = fopen(logName, "w");
      else {
        status = STATUS_LOG_OPEN_FAILED;
        continue;
      }

      // If log file is not opened then signal an error
      if (gLogFilePtr == NULL) {
        status = STATUS_LOG_OPEN_FAILED;
        continue;
      }

      // Store log task name
      sprintf(gLogTaskName, "%s", taskName);

    } while (0); // End of main do-loop

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Finish task logging
 *
 * @param[in] status Error status.
 ******************************************************************************/
Status LogClose(Status status) {

    // Declare (and initialise) variables

    // Main do-loop to fall through in case of an error
    do {

      // Close log file
      if (gLogFilePtr != NULL)
        fclose(gLogFilePtr);
      else {
        status = STATUS_LOG_CLOSE_FAILED;
        continue;
      }

    } while (0); // End of main do-loop

    // Return status
    return status;

}


/**************************************************************************//**
 * @brief Log message
 *
 * @param[in] msgType Message type.
 * @param[in] msgFormat Pointer to message format.
 ******************************************************************************/
Status Log(MessageType msgType, const char *msgFormat, ...) {

    // Declare (and initialise) variables
    Status    status = STATUS_OK;
    va_list   vl;
    time_t    now;
    struct tm timeStruct;
    char      type[100];

    // Main do-loop to fall through in case of an error
    do {

      // If no log file has been opened then open one now
      if (gLogFilePtr == NULL) {
        status = LogInit(DEFAULT_LOG_FILENAME, DEFAULT_TASK_NAME, status);
        if (status != STATUS_OK)
          continue;
      }

      // Set message type
      switch (msgType) {
      case Error_0:
        sprintf(type, "Error_0");
        break;
      case Error_1:
        sprintf(type, "Error_1");
        break;
      case Error_2:
        sprintf(type, "Error_2");
        break;
      case Error_3:
        sprintf(type, "Error_3");
        break;
      case Warning_0:
        sprintf(type, "Warn_0 ");
        break;
      case Warning_1:
        sprintf(type, "Warn_1 ");
        break;
      case Warning_2:
        sprintf(type, "Warn_2 ");
        break;
      case Warning_3:
        sprintf(type, "Warn_3 ");
        break;
      case Alert_0:
        sprintf(type, "Alert_0");
        break;
      case Alert_1:
        sprintf(type, "Alert_1");
        break;
      case Alert_2:
        sprintf(type, "Alert_2");
        break;
      case Alert_3:
        sprintf(type, "Alert_3");
        break;
      case Log_0:
        sprintf(type, "Log_0  ");
        break;
      case Log_1:
        sprintf(type, "Log_1  ");
        break;
      case Log_2:
        sprintf(type, "Log_2  ");
        break;
      case Log_3:
        sprintf(type, "Log_3  ");
        break;
      default:
        sprintf(type, "Unknown");
        break;
      }

      // Get time
      now = time(NULL);
      #ifdef HAVE_GMTIME_R
      gmtime_r(&now, &timeStruct);
      #else
      memcpy(&timeStruct, gmtime(&now), sizeof(struct tm));
      #endif

      // Write message type, time and task name to log file
      if (fprintf(gLogFilePtr, "%s %04d-%02d-%02dT%02d:%02d:%02d %s: ",
                  type,
                  timeStruct.tm_year + 1900,
                  timeStruct.tm_mon + 1,
                  timeStruct.tm_mday,
                  timeStruct.tm_hour,
                  timeStruct.tm_min,
                  timeStruct.tm_sec,
                  gLogTaskName) < 0) {
        status = STATUS_LOG_WRITE_FAILED;
        continue;
      }

      // Write message to log file
      va_start(vl, msgFormat);
      if (vfprintf(gLogFilePtr, msgFormat, vl) < 0) {
        status = STATUS_LOG_WRITE_FAILED;
        continue;
      }
      va_end(vl);

      // Write <CR> to log file
      if (fprintf(gLogFilePtr, "\n") != 1) {
        status = STATUS_LOG_WRITE_FAILED;
        continue;
      }

    } while (0); // End of main do-loop

    // Return status
    return status;

}

/* Namespace ends ___________________________________________________________ */
}
