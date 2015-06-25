HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_CXXTASK              = gtobssim

HD_CXXTASK_SRC_cxx      = obsSim.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_OBSERVATIONSIM_LIBS} \
                          ${GT_FITSGEN_LIBS} \
                          -levtUtils \
                          ${GT_DATASUBSELECTOR_LIBS} \
                          ${GT_EVTBIN_LIBS} \
                          ${GT_SOURCE_LIBS} \
                          -lirfLoader \
                          -lirfInterface \
                          -lirfUtil \
                          -llatResponse \
                          -ldc1aResponse \
                          -ltestResponse \
                          ${GT_XML_LIBS} \
                          ${GT_EMBED_PYTHON_LIBS} \
                          -lTreePlayer \
                          ${PYTHON_LIB} \
                          ${GT_STD_CXXLIBS}

HD_INSTALL_TASKS        = ${HD_CXXTASK}

HD_INSTALL_HELP         = gtobssim.txt

include ${HD_STD_MAKEFILE}
