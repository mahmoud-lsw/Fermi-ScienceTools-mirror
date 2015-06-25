HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_TEST_TASKS           = test_EarthPhenom

HD_TEST_CXXTASK         = ${HD_TEST_TASKS}

HD_TEST_CXXTASK_SRC_cxx	= test_main.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          -lEarthPhenom \
                          ${GT_FLUX_LIBS} \
                          ${GT_XML_LIBS} \
                          ${GT_STD_CXXLIBS}

include ${HD_STD_MAKEFILE}