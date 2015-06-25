/** \file st_app_main.cxx
    \brief Standard main for ballistic Science Tools applications.
    \author James Peachey, HEASARC
*/

// Add floating point exception traps if requested.
#ifdef BUILD_WITHOUT_ROOT
#include <Python.h>
#endif

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <iostream>
#include <stdexcept>
#include <typeinfo>

#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "facilities/commonUtilities.h"
#ifdef WIN32
#include "facilities/AssertDialogOverride.h"
#endif

int main(int argc, char ** argv) {
#ifdef BUILD_WITHOUT_ROOT
	Py_Initialize();
	PySys_SetArgv(argc, argv);
#endif

#ifdef _DEBUG
   _CrtSetReportHook( AssertDialogOverride );
   _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
   _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
   _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
#endif

  int status = 0;
  st_app::StApp * this_st_app = 0;

  bool debug = false;
  bool gui_mode = false;

  try {
#ifdef TRAP_FPE
    // Add floating point exception traps.
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
    //Set up the necessary environment variables to be used by applications
    facilities::commonUtilities::setupEnvironment();

    // Process command line arguments. This will throw if the real application code
    // cannot/should not start.
    st_app::StApp::processCommandLine(argc, argv);

    // Read standard parameters to set up chatter, debug mode, etc.
    st_app::IStAppFactory::instance().configureApp();

    // Determine whether debugging was selected, either at runtime or compile time.
    debug = st_app::IStAppFactory::instance().getDebugMode();

    // Determine whether to run the GUI.
    gui_mode = st_app::IStAppFactory::instance().getGuiMode();

    if (!debug) {
      // Not in debug mode, so create and run the application here, inside the top-level try block.
      try {
        // Try using the singleton StAppFactory to create the application:
        this_st_app = st_app::IStAppFactory::instance().createApp();

        // Set application name.
        this_st_app->setName(st_app::IStAppFactory::instance().getAppName());
      } catch(const std::logic_error &) {
        throw std::logic_error("Failed to get an StAppFactory singleton: client must define one");
      }

      // Check debug status again in case the application constructor changed it, and defer running if
      // debug mode is currently selected.
      if (0 != this_st_app && !debug) {
        // Display startup banner.
        this_st_app->banner();

        // Run the application, either in GUI mode or not.
        if (gui_mode) this_st_app->runGui();
        else this_st_app->run();
      }
    }

  } catch (const std::exception & x) {
    // Return a non-zero exit code:
    status = 1;

    // Report the type of the exception if possible, using typeid; typeid can throw so be careful:
    const char * type_name = "std::exception";
    try {
      type_name = typeid(x).name();
    } catch (...) {
      // Ignore problems with typeid.
    }
    std::cerr << "Caught " << type_name << " at the top level: " << x.what() << std::endl;
  }

  if (debug) {
    // In debug mode, so create and run the application here, outside of any try-catch,
    // if it was not already created,
    // so that the debugger is more useful for finding exceptions at their source.
    // Use the singleton StAppFactory to create the application:
    if (0 == this_st_app) {
      this_st_app = st_app::IStAppFactory::instance().createApp();

      // Set application name.
      this_st_app->setName(st_app::IStAppFactory::instance().getAppName());
    }

    if (0 != this_st_app) {
      // Display startup banner.
      this_st_app->banner();

      // Run the application, either in GUI mode or not.
      if (gui_mode) this_st_app->runGui();
      else this_st_app->run();
    }
  }

  // Clean up:
  delete this_st_app;

  return status;
}
