HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} ${GT_FACILITIES_LIBS}

HD_TEST_TASKS           = test_env test_time testUtil

HD_CLEAN		= test_env${BSUF} test_time${BSUF} testUtil${BSUF}

include ${HD_STD_MAKEFILE}

test:
	${HD_MAKE} test_env HD_TEST_CXXTASK=test_env HD_TEST_CXXTASK_SRC_cxx=test_env.cxx
	${HD_MAKE} test_time HD_TEST_CXXTASK=test_time HD_TEST_CXXTASK_SRC_cxx=test_time.cxx
	${HD_MAKE} testUtil HD_TEST_CXXTASK=testUtil HD_TEST_CXXTASK_SRC_cxx=testUtil.cxx
