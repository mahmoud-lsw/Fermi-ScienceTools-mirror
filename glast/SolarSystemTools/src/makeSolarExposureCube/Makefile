HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_CXXTASK              = gtltcubesun

HD_CXXTASK_SRC_cxx      = makeSolarExposureCube.cxx

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

HD_INSTALL_HELP         = gtltcubesun.txt

include ${HD_STD_MAKEFILE}
