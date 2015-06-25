HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = _pyIrfLoader

HD_LIBRARY_SRC_cxx      = pyIrfLoader_wrap.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS} ${PYTHON_INC} -fpermissive

HD_ADD_SHLIB_LIBS       = yes

SHLIB_LD_CXX            = ${PYLIKE_SHLIB_LD_CXX}

SHLIB_SUFFIX            = ${PYLIKE_SHLIB_SUFFIX}

HD_SHLIB_LIBS           = ${HD_LFLAGS} \
                          -l${CLHEP} \
                          -ldc1aResponse \
                          -lirfLoader \
                          -lirfInterface \
                          -lastro \
                          -lfacilities \
                          -ltip \
                          -lst_facilities \
                          -lf2c \
                          -llatResponse \
                          -ltestResponse \
                          -lst_stream \
                          -lirfUtil \
                          ${ROOT_GUI_LIBS} \
                          ${HD_STD_LIBS} \
                          ${SYSLIBS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_INSTALL_PYTHON_LIBS  = pyIrfLoader.py

HD_DISTCLEAN            = pyIrfLoader_wrap.cxx pyIrfLoader.py

default: custom-lib

all: default publish

include ${HD_STD_MAKEFILE}

custom-lib:
	@if [ "x${BUILD_WITHOUT_ROOT}" = xno ]; then \
		${HD_MAKE} ${HD_LIBRARY_ROOT}; \
	else \
		${HD_MAKE} ${HD_LIBRARY_ROOT} CERN_ROOT_CXXLIBS=""; \
	fi
	${HD_MAKE} ${HD_LIBRARY_ROOT}

pyIrfLoader_wrap.cxx: pyIrfLoader.i
	@PATH=${EXT_BLD_EXEC_PFX}/bin:${EXT_BLD_EXEC_PFX}/lib:${HD_TOP_BIN}:${HD_TOP_LIB}:${PATH}; \
	LD_LIBRARY_PATH=${EXT_BLD_EXEC_PFX}/lib:${HD_TOP_LIB}:${LD_LIBRARY_PATH}; \
	DYLD_LIBRARY_PATH=${EXT_BLD_EXEC_PFX}/lib:${HD_TOP_LIB}:${DYLD_LIBRARY_PATH}; \
	ROOTSYS=${shell ${HD_BUILD_DIR}/hd_scanenv ${HD_SRC}/../external/BUILD_DIR/hmakerc HD_BLD_EXEC_PFX}; \
	export PATH; \
	export LD_LIBRARY_PATH; \
	export DYLD_LIBRARY_PATH; \
	export ROOTSYS; \
	SWIG_VERSION=${shell echo $${SWIG} | sed s/swig//}; \
	export SWIG_VERSION; \
	swig -c++ -python \
	-I${HD_BLD_EXEC_PFX}/include -I${EXT_BLD_EXEC_PFX}/include \
	-I${EXT_BLD_EXEC_PFX}/share/swig/$${SWIG_VERSION}/python \
	-I${EXT_BLD_EXEC_PFX}/share/swig/$${SWIG_VERSION} -I.. \
	-interface lib_pyIrfLoader pyIrfLoader.i
