HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = st_graph

CERN_ROOT_SRC_cxx       = RootEngine.cxx \
                          RootFrame.cxx \
                          RootPlot.cxx \
                          RootPlotFrame.cxx \
                          RootTabFolder.cxx \
                          STGMainFrame.cxx \
                          STGRootCint.cxx
                          
MPL_SRC_cxx				= EmbedPython.cxx \
						  MPLEngine.cxx \
						  MPLFrame.cxx \
						  MPLPlotFrame.cxx \
						  MPLPlot.cxx \
						  MPLTabFolder.cxx

HD_LIBRARY_SRC_cxx      = Axis.cxx \
                          Engine.cxx \
                          IPlot.cxx \
                          StGui.cxx \
                          ${MPL_SRC_cxx} \
                          ${CERN_ROOT_SRC_cxx}

HD_CXXFLAGS             = ${PYTHON_INC} ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_INSTALL_PYTHON_LIBS  = STTopLevel.py \
                          Lego.py

HD_CLEAN                = STGRootCint.cxx STGRootCint.h

default: custom-lib default-subdirs

all: default publish all-subdirs

include ${HD_STD_MAKEFILE}

custom-lib:
	@if [ "x${BUILD_WITHOUT_ROOT}" = xno ]; then \
		${HD_MAKE} ${HD_LIBRARY_ROOT} MPL_SRC_cxx=""; \
	else \
		${HD_MAKE} ${HD_LIBRARY_ROOT} CERN_ROOT_SRC_cxx=""; \
	fi

STGRootCint.cxx: ../st_graph/RootFrame.h ../st_graph/LinkDef.h
	@PATH=${EXT_BLD_EXEC_PFX}/bin:${EXT_BLD_EXEC_PFX}/lib:${HD_TOP_BIN}:${HD_TOP_LIB}:${PATH}; \
	LD_LIBRARY_PATH=${EXT_BLD_EXEC_PFX}/lib:${HD_TOP_LIB}:${LD_LIBRARY_PATH}; \
	DYLD_LIBRARY_PATH=${EXT_BLD_EXEC_PFX}/lib:${HD_TOP_LIB}:${DYLD_LIBRARY_PATH}; \
	ROOTSYS=${shell ${HD_BUILD_DIR}/hd_scanenv ${HD_SRC}/../external/BUILD_DIR/hmakerc HD_BLD_EXEC_PFX}; \
	export PATH; \
	export LD_LIBRARY_PATH; \
	export DYLD_LIBRARY_PATH; \
	export ROOTSYS; \
	rootcint -f $@ -c -I. -I.. -I$$ROOTSYS/include \
	../st_graph/RootFrame.h ../st_graph/LinkDef.h
