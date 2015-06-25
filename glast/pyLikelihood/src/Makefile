HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = _pyLikelihood

HD_LIBRARY_SRC_cxx      = pyLikelihood_wrap.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS} ${PYTHON_INC} -fpermissive

HD_ADD_SHLIB_LIBS       = yes

CERN_ROOT_CXXLIBS       = -lCore -lCint -lRIO -lNet -lTree -lRint -lMatrix \
                          -lHist -lGraf -lGpad -lPostscript -lMatrix \
                          -lPhysics -lGui -lTreePlayer -lProof -lGX11 \
                          -lGX11TTF -lHist -lGraf -lGraf3d -lGpad -lThread \
                          -lHistPainter

SHLIB_LD_CXX            = ${PYLIKE_SHLIB_LD_CXX}

SHLIB_SUFFIX            = ${PYLIKE_SHLIB_SUFFIX}

HD_SHLIB_LIBS           = ${HD_LFLAGS} \
                          -lLikelihood -lastro -ldataSubselector -lirfLoader -ldc1aResponse \
                          -leblAtten -levtbin -lfacilities -lhealpix -lirfInterface \
                          -lirfUtil \
                          -llatResponse -lmap_tools -loptimizers -lst_app \
                          -lst_facilities -lst_graph -lst_stream \
                          -ltestResponse -ltip -lxmlBase \
                          -l${CLHEP} -l${CPPUNIT} -lfftw3 ${PYTHON_LIB} \
                          ${ROOT_GUI_LIBS} -lMinuit2 -l${XERCES} \
                          ${CERN_ROOT_CXXLIBS} \
                          ${HD_STD_LIBS} ${SYSLIBS} -lf2c

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_INSTALL_PYTHON_LIBS  = pyLikelihood.py

HD_DISTCLEAN            = pyLikelihood_wrap.cxx pyLikelihood.py

default: custom-lib

all: default publish

include ${HD_STD_MAKEFILE}

custom-lib:
	@if [ "x${BUILD_WITHOUT_ROOT}" = xno ]; then \
		${HD_MAKE} ${HD_LIBRARY_ROOT}; \
	else \
		${HD_MAKE} ${HD_LIBRARY_ROOT} CERN_ROOT_CXXLIBS=""; \
	fi

pyLikelihood_wrap.cxx: pyLikelihood.i
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
	-I${HD_BLD_EXEC_PFX}/include \
	-I${EXT_BLD_EXEC_PFX}/share/swig/$${SWIG_VERSION}/python \
	-I${EXT_BLD_EXEC_PFX}/share/swig/$${SWIG_VERSION} -I.. \
	-interface lib_pyLikelihood pyLikelihood.i
