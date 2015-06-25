HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = GRB

HD_LIBRARY_SRC_cxx      = GRB/GRBConstants.cxx \
                          GRB/GRBengine.cxx \
                          GRB/GRBmanager.cxx \
                          GRB/GRBShell.cxx \
                          GRB/GRBShock.cxx \
                          GRB/GRBSim.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_CLEAN                = ${HD_LIBRARY_SRC_cxx:.cxx=.o}

include ${HD_STD_MAKEFILE}
