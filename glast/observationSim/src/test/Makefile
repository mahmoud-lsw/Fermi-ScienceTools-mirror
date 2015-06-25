HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_TEST_TASKS           = test_observationSim

HD_TEST_CXXTASK         = ${HD_TEST_TASKS}

HD_TEST_CXXTASK_SRC_cxx	= main.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_OBSERVATIONSIM_LIBS} \
                          ${GT_FITSGEN_LIBS} \
                          -levtUtils \
                          -lTreePlayer \
                          ${GT_DATASUBSELECTOR_LIBS} \
                          ${GT_EVTBIN_LIBS} \
                          ${GT_SOURCE_LIBS} \
                          ${GT_FLUX_LIBS} \
                          -lirfLoader \
                          -lirfInterface \
                          -lirfUtil \
                          -llatResponse \
                          -ldc1aResponse \
                          -ltestResponse \
                          ${GT_XML_LIBS} \
                          ${GT_EMBED_PYTHON_LIBS} \
                          ${PYTHON_LIB} \
                          ${GT_STD_CXXLIBS}

include ${HD_STD_MAKEFILE}
