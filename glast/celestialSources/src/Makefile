HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = celestialSources

CERN_ROOT_SRC_cxx       = TRandom4.cxx

HD_LIBRARY_SRC_cxx      = SpectrumFactoryLoader.cxx ${CERN_ROOT_SRC_cxx}

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

default: custom-lib default-subdirs

all: default publish all-subdirs

include ${HD_STD_MAKEFILE}

custom-lib:
	@if [ "x${BUILD_WITHOUT_ROOT}" = xno ]; then \
		${HD_MAKE} ${HD_LIBRARY_ROOT}; \
	else \
		${HD_MAKE} ${HD_LIBRARY_ROOT} CERN_ROOT_SRC_cxx=""; \
	fi
