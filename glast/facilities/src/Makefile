HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = facilities

HD_LIBRARY_SRC_cxx      = Clock.cxx \
                          ScheduledEvent.cxx \
                          Scheduler.cxx \
                          Timestamp.cxx \
                          Util.cxx \
                          XGetopt.cxx \
                          commonUtilities.cxx

# The -DHEADAS is needed to ensure some HEADAS-conditional code is
# compiled in ccommonUtilities.cxx.
HD_CXXFLAGS             = ${HD_STD_CXXFLAGS} -DHEADAS

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

include ${HD_STD_MAKEFILE}
