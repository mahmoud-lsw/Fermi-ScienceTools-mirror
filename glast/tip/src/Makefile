HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = tip

CERN_ROOT_SRC_cxx       = RootTable.cxx

HD_LIBRARY_SRC_cxx      = FileSummary.cxx \
                          FitsFileManager.cxx \
                          FitsHeader.cxx \
                          FitsPrimProps.cxx \
                          FitsTable.cxx \
                          FitsTipFile.cxx \
                          Header.cxx \
                          IFileSvc.cxx \
                          KeyRecord.cxx \
                          LinearInterp.cxx \
                          TipException.cxx \
                          TipFile.cxx \
                          ${CERN_ROOT_SRC_cxx} 

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

# <HACK>
# Need to install this local header to get Likelihood to compile.
HD_INC_SUBDIR           = src
HD_INSTALL_HEADERS      = FitsTipFile.h
# </HACK>

HD_INSTALL_LIBRARIES	= ${HD_LIBRARY_ROOT}

default: custom-lib default-subdirs

all: default publish all-subdirs

include ${HD_STD_MAKEFILE}

custom-lib:
	@if [ "x${BUILD_WITHOUT_ROOT}" = xno ]; then \
		${HD_MAKE} ${HD_LIBRARY_ROOT}; \
	else \
		${HD_MAKE} ${HD_LIBRARY_ROOT} CERN_ROOT_SRC_cxx=""; \
	fi
