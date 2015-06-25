HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       = 

HD_CXXTASK              = gtexposure

HD_CXXTASK_SRC_cxx      = gtexposure.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_STD_CXXLIBS} \
                          ${GT_DATASUBSELECTOR_LIBS} \
                          ${GT_OPT_LIBS} \
                          ${GT_XML_LIBS} \
                          ${GT_EVTBIN_LIBS} \
                          ${GT_PYEXPOSURE_LIBS} \
                          ${GT_LIKELIHOOD_LIBS} \
                          ${GT_HEALPIX_LIBS} \
                          ${GT_MAP_TOOLS_LIBS} \
                          ${GT_SOURCE_LIBS} \
                          ${GT_FFT_LIBS} \
                          -lirfLoader \
                          -lirfInterface \
                          -lirfUtil \
                          -llatResponse \
                          -ldc1aResponse \
                          -ltestResponse

HD_INSTALL_HELP         = gtexposure.txt

HD_INSTALL_TASKS        = ${HD_CXXTASK}

include ${HD_STD_MAKEFILE}
