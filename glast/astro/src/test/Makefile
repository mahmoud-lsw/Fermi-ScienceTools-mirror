HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_TEST_TASKS           = test_astro

HD_TEST_CXXTASK         = ${HD_TEST_TASKS}

HD_TEST_CXXTASK_SRC_cxx	= test.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_ASTRO_LIBS} \
                          ${GT_TIP_LIBS} \
                          ${GT_FACILITIES_LIBS} \
                          ${GT_CLHEP_LIBS} \
                          ${ROOT_GUI_LIBS} \
                          ${HD_STD_LIBS} \
                          ${SYSLIBS}

include ${HD_STD_MAKEFILE}

#Do not change this.  It is needed to work around the behavior
#of make when using the test target and a file named test.cxx
test:
	${HD_MAKE} ${HD_TEST_CXXTASK}
