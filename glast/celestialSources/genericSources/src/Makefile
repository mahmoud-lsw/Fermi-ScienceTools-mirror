HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = genericSources

HD_LIBRARY_SRC_c        = dgaus8.c

CERN_ROOT_SRC_cxx       = TF1Map.cxx TF1Spectrum.cxx

HD_LIBRARY_SRC_cxx      = FileSpectrum.cxx \
                          FileSpectrumMap.cxx \
                          FitsImage.cxx \
                          FitsTransient.cxx \
                          GaussianSource.cxx \
                          GaussianSpectrum.cxx \
                          Isotropic.cxx \
                          IsotropicFileSpectrum.cxx \
                          MapCube.cxx \
                          MapSource.cxx \
                          PeriodicSource.cxx \
                          Pulsar.cxx \
                          RadialSource.cxx \
                          SimpleTransient.cxx \
                          SourcePopulation.cxx \
                          SpectralTransient.cxx \
                          TransientTemplate.cxx \
                          Util.cxx \
                          ${CERN_ROOT_SRC_cxx}

HD_CFLAGS               = ${HD_STD_CFLAGS}

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

default: custom-lib default-subdirs

all: default publish all-subdirs

include ${HD_STD_MAKEFILE}

custom-lib:
	@if [ "x${BUILD_WITHOUT_ROOT}" = xno ]; then \
		${HD_MAKE} ${HD_LIBRARY_ROOT}; \
	else \
		${HD_MAKE} ${HD_LIBRARY_ROOT} CERN_ROOT_SRC_cxx=""; \
	fi

Util.o: Util.cxx
	${CXX} -c  -o ${<D}/${@F} ${HD_CXXFLAGS} -O0 $<