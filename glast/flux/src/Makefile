HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = flux

HD_LIBRARY_SRC_cxx      = CompositeSource.cxx \
                          Earth.cxx \
                          EventSource.cxx \
                          FileSource.cxx \
                          Flux.cxx \
                          FluxMgr.cxx \
                          FluxSource.cxx \
                          SimpleSpectrum.cxx \
                          SourceDirection.cxx \
                          Spectrum.cxx \
                          SpectrumFactoryTable.cxx \
                          SurfaceMuons.cxx \
                          TimeCandle.cxx \
                          VdgGamma.cxx \
                          rootplot/rootAngleHist.cxx \
                          rootplot/rootEnergyHist.cxx \
                          rootplot/rootHist.cxx \
                          rootplot/rootplot.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_CLEAN                = ${HD_LIBRARY_SRC_cxx:.cxx=.o}

include ${HD_STD_MAKEFILE}
