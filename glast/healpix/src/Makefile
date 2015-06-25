HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = healpix

HD_LIBRARY_SRC_cxx      = AlmOp.cxx \
                          CosineBinner.cxx \
                          HealpixArrayIO.cxx \
                          Healpix.cxx \
                          HealPixel.cxx \
                          HealpixMap.cxx \
                          Map.cxx

HD_LIBRARY_SRC_cc       = base/alm_filter_tools.cc \
                          base/alm_map_tools.cc \
                          base/alm_powspec_tools.cc \
                          base/cxxutils.cc \
                          base/fitshandle.cc \
                          base/healpix_base.cc \
                          base/healpix_map.cc \
                          base/healpix_map_fitsio.cc \
                          base/powspec.cc \
                          base/rotmatrix.cc

HD_LIBRARY_SRC_c        = base/bluestein.c \
                          base/fftpack.c \
                          base/ls_fft.c

HD_SUBDIRS              = base

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS} -fpermissive

HD_CFLAGS               = ${HD_STD_CFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_CLEAN                = ${HD_LIBRARY_SRC_cxx:.cxx=.o} \
                          ${HD_LIBRARY_SRC_cx:.cc=.o} \
                          ${HD_LIBRARY_SRC_c:.c=.o}

include ${HD_STD_MAKEFILE}
