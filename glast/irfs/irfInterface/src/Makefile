HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = irfInterface

HD_LIBRARY_SRC_cxx      = IEdisp.cxx \
                          IEfficiencyFactor.cxx \
                          IPsf.cxx \
                          IrfRegistry.cxx \
                          IrfsFactory.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_EDISP_CXXFLAGS       = ${HD_CXXFLAGS} -O0

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

include ${HD_STD_MAKEFILE}

IEdisp.o: IEdisp.cxx
	${CXX} -c  -o $@ ${HD_EDISP_CXXFLAGS} $?
