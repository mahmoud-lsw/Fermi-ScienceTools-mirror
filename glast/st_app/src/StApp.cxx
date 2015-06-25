/** \file StApp.cxx
    \brief Base class for Science Tools applications.
    \author James Peachey, HEASARC
*/
#include <cctype>
#include <stdexcept>
#include <string>

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppGui.h"
#include "st_app/StGui.h"

#include "st_graph/Engine.h"

#include "st_stream/StreamFormatter.h"

namespace st_app {

  // Static member definitions:
  int StApp::s_argc = 0;
  char ** StApp::s_argv = 0;

  // Process (for now simply store) command line arguments:
  void StApp::processCommandLine(int argc, char ** argv) {
    s_argc = argc;
    s_argv = argv;
  }

  // Accessor for number of command line arguments.
  int StApp::getArgc() { return s_argc; }

  // Accessor for command line arguments.
  char ** StApp::getArgv() { return s_argv; }

  // Construct application object:
  StApp::StApp(): m_name("Unknown application"), m_version(), m_par_group(0), m_gui(0), m_main_frame(0), m_plot_frame(0) {
    setVersion("");
  }

  // Destruct application object:
  StApp::~StApp() throw() { delete m_main_frame; delete m_gui; delete m_par_group; }

  void StApp::runGui() {
    // TODO Clean this up when old GUI removed. For now it's a rats' nest.
    // If no GUI already instantiated, create the old GUI.
    if (0 == m_gui) {
      StEventReceiver * gui = new StEventReceiver(st_graph::Engine::instance(), getParGroup(), this);
      m_gui = gui;
      gui->run();
    } else {
      // See if this is the new GUI (set by client) and run that if it is.
      st_graph::StGui * gui = dynamic_cast<st_graph::StGui *>(m_gui);
      if (0 != gui) {
        gui->run();
      } else {
        // See if this is the old GUI (set by client) and run that if it is.
        StEventReceiver * gui = dynamic_cast<StEventReceiver *>(m_gui);
        if (0 != gui) {
          gui->run();
        } else {
          // Unknown type of GUI (set by client) so just run the event loop.
          st_graph::Engine::instance().run();
        }
      }
    }
  }

  // Return an object which provides the most typical kinds of parameter access.
  AppParGroup & StApp::getParGroup() {
    return getParGroup(m_name);
  }

  // Return an object which provides the most typical kinds of parameter access.
  AppParGroup & StApp::getParGroup(const std::string & app_name) {
    // Create if necessary:
    if (0 == m_par_group) m_par_group = new AppParGroup(app_name);
    return *m_par_group;
  }

  void StApp::banner() const {
    st_stream::StreamFormatter sf("StApp", "banner", 2);
    sf.info(1) << "This is " << m_name << " version " << m_version << std::endl;
  }

  const std::string & StApp::getName() const { return m_name; }

  const std::string & StApp::getVersion() const { return m_version; }

  void StApp::setName(const std::string & name) { m_name = name; }

  void StApp::setVersion(const std::string & version) {
    static const std::string cvs_prefix("$Name:");

    std::string::size_type start = version.find(cvs_prefix);

    std::string default_version = "N/A";

    // If cvs was used to assign the version automatically, chop out the prefix and suffix ($) it uses.
    if (std::string::npos != start) {
      // Change default version to reflect the fact that a CVS prefix was found.
      default_version = "HEAD";

      start += cvs_prefix.size();

      // Skip whitespace after prefix.
      for (std::string::const_iterator itor = version.begin() + start; itor != version.end(); ++itor, ++start) {
        if (!isspace(*itor)) break;
      }

      // Stop when whitespace or trailing $ is reached.
      std::string::size_type stop = start;
      for (std::string::const_iterator itor = version.begin() + stop; itor != version.end(); ++itor, ++stop) {
        if (isspace(*itor) || '$' == *itor) break;
      }

      // String between start and stop is the version.
      m_version = version.substr(start, stop - start);

    } else {
      m_version = version;
    }

    // See if version contains any non-whitespace.
    std::string::const_iterator vers_itor;
    for (vers_itor = m_version.begin(); vers_itor != m_version.end() && (0 != isspace(*vers_itor)); ++vers_itor) {}

    // If no non-whitespace was found, call this version HEAD.
    if (vers_itor == m_version.end()) m_version = default_version;
  }

  st_graph::IFrame * StApp::getPlotFrame(const std::string & title) {
    st_graph::IFrame * plot_frame(0);
    st_graph::StGui * st_gui = dynamic_cast<st_graph::StGui *>(m_gui);

    if (0 == st_gui) {
      st_graph::Engine & engine(st_graph::Engine::instance());
      if (0 == m_main_frame) m_main_frame = engine.createMainFrame(0, 638, 319, "StApp plot window");
      if (0 == m_plot_frame) m_plot_frame = engine.createPlotFrame(m_main_frame, title, 638, 319);
      plot_frame = m_plot_frame;
    } else {
      // In a GUI context use the Gui's plot window, creating it as needed.
      plot_frame = st_gui->getPlotFrame();
    }

    return plot_frame;
  }
}
