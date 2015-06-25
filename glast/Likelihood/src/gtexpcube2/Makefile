HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_CXXTASK              = gtexpcube2

HD_CXXTASK_SRC_cxx      = gtexpcube2.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          -lLikelihood \
                          -loptimizers \
                          -levtbin \
                          -ldataSubselector \
                          -lhealpix \
                          -lmap_tools \
                          -lxmlBase \
                          -lirfLoader \
                          -lirfInterface \
                          -leblAtten \
                          -llatResponse \
                          -ldc1aResponse \
                          -ltestResponse \
                          -lirfUtil \
                          -lxerces-c \
                          -lfftw3 \
                          -lMinuit2 \
                          ${GT_STD_CXXLIBS}

HD_INSTALL_TASKS        = ${HD_CXXTASK}

HD_INSTALL_HELP         = gtexpcube2.txt

include ${HD_STD_MAKEFILE}
