HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_TEST_TASKS           = test_fitsGen

HD_TEST_CXXTASK         = ${HD_TEST_TASKS}

HD_TEST_CXXTASK_SRC_cxx	= test.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS} -I../..

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_FITSGEN_LIBS} \
                          -levtUtils \
                          -lxmlBase \
                          -lxerces-c \
                          -lTreePlayer \
                          ${GT_EVTBIN_LIBS} \
                          ${GT_EMBED_PYTHON_LIBS} \
			  ${PYTHON_LIB} \
                          ${GT_STD_CXXLIBS}

include ${HD_STD_MAKEFILE}

#Do not change this.  It is needed to work around the behavior
#of make when using the test target and a file named test.cxx
test:
	${HD_MAKE} ${HD_TEST_CXXTASK}
