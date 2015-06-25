HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = irfUtil

HD_LIBRARY_SRC_c        = HDgtcalf.c HDgtcalf_internal.c HDgtcalf_standalone.c

CERN_ROOT_SRC_cxx       = RootTable.cxx

HD_LIBRARY_SRC_cxx      = HdCaldb.cxx Util.cxx IrfHdus.cxx ${CERN_ROOT_SRC_cxx}

HD_CFLAGS               = ${HD_STD_CFLAGS}

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

default: custom-lib

all: default publish

include ${HD_STD_MAKEFILE}

custom-lib:
	@if [ "x${BUILD_WITHOUT_ROOT}" = xno ]; then \
		${HD_MAKE} ${HD_LIBRARY_ROOT}; \
	else \
		${HD_MAKE} ${HD_LIBRARY_ROOT} CERN_ROOT_SRC_cxx=""; \
	fi
