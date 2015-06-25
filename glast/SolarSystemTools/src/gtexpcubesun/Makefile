HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_CXXTASK              = gtexpcubesun

HD_CXXTASK_SRC_cxx      = gtexpcubesun.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          -lLikelihood \
                          -loptimizers \
                          -levtbin \
                          -lhealpix \
                          -ldataSubselector \
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
                          -lSolarSystemTools \
                          ${GT_STD_CXXLIBS}

HD_INSTALL_TASKS        = ${HD_CXXTASK}

include ${HD_STD_MAKEFILE}
