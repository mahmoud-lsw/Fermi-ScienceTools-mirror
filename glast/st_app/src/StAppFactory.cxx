/** \file StAppFactory.cxx
    \brief Factory class for Science Tools application objects derived from StApp.
    \author James Peachey, HEASARC
*/

#include <stdexcept>

#include "hoops/hoops_exception.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_stream/Stream.h"
#include "st_stream/st_stream.h"

namespace st_app {

  // Static member definitions:
  IStAppFactory * IStAppFactory::s_factory = 0;

  // Singleton access.
  IStAppFactory & IStAppFactory::instance() {
    if (0 == s_factory) throw std::logic_error("IStAppFactory::instance() cannot find a factory");
    return *s_factory;
  }

  // Construct IStAppFactory, setting singleton if it's not already defined.
  IStAppFactory::IStAppFactory(): m_app_name(), m_gui_mode(false) { if (0 == s_factory) s_factory = this; }

  // Construct IStAppFactory, setting singleton if it's not already defined.
  IStAppFactory::IStAppFactory(const std::string & app_name): m_app_name(app_name), m_gui_mode(false) {
    if (0 == s_factory) s_factory = this;
  }

  // Destruct IStAppFactory, unsetting singleton if this is it.
  IStAppFactory::~IStAppFactory() throw() { if (this == s_factory) s_factory = 0; }

  void IStAppFactory::configureApp() {
    // Initialize standard st_stream objects.
    st_stream::OStream::initStdStreams();

    // If application has a name, use it to set up various standard global properties.
    if (!m_app_name.empty()) {
      // Pass the name to st_stream, to use in its prefixes.
      st_stream::SetExecName(m_app_name);
      try {
        // Get parameter object, ignoring problems.
        st_app::AppParGroup pars(m_app_name);
        try {
          // Try to get chatter parameter, and use it to set the maximum tool chatter. Ignore problems.
          int chat = pars["chatter"];
          setMaximumChatter(chat);
        } catch (hoops::Hexception &) {}
        try {
#ifdef ST_APP_DEBUG
          // Compile-time debugging trumps debug parameter.
          setDebugMode(true);
#else
          // Try to get debug parameter, and use it to enable/disable debug behavior. Ignore problems.
          bool debug_mode = pars["debug"];
          setDebugMode(debug_mode);
#endif
        } catch (hoops::Hexception &) {}
        try {
          // Try to get gui parameter, and use it to enable/disable Gui. Ignore problems.
          m_gui_mode = pars["gui"];
        } catch (hoops::Hexception &) {}
      } catch (hoops::Hexception &) {}
    }
  }

  int IStAppFactory::getMaximumChatter() const { return st_stream::GetMaximumChatter(); }

  void IStAppFactory::setMaximumChatter(int maximum_chatter) { st_stream::SetMaximumChatter(maximum_chatter); }

  bool IStAppFactory::getDebugMode() const { return st_stream::GetDebugMode(); }

  void IStAppFactory::setDebugMode(bool debug_mode) { st_stream::SetDebugMode(debug_mode); }

  bool IStAppFactory::getGuiMode() const { return m_gui_mode; }

  void IStAppFactory::setGuiMode(bool gui_mode) { m_gui_mode = gui_mode; }

  const std::string & IStAppFactory::getAppName() const { return m_app_name; }

}
