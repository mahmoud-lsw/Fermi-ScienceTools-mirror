HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = astro

HD_LIBRARY_SRC_c        = wcslib/cel.c \
                          wcslib/fitshdr.c \
                          wcslib/getwcstab.c \
                          wcslib/lin.c \
                          wcslib/log.c \
                          wcslib/prj.c \
                          wcslib/spc.c \
                          wcslib/sph.c \
                          wcslib/spx.c \
                          wcslib/tab.c \
                          wcslib/wcs.c \
                          wcslib/wcsfix.c \
                          wcslib/wcshdr.c \
                          wcslib/wcspih.c \
                          wcslib/wcstrig.c \
                          wcslib/wcsulex.c \
                          wcslib/wcsunits.c \
                          wcslib/wcsutil.c \
                          wcslib/wcsutrn.c

HD_LIBRARY_SRC_cxx      = EarthCoordinate.cxx \
                          EarthOrbit.cxx \
                          Geomag.cxx \
                          GPS.cxx \
                          HTM.cxx \
                          IGRField.cxx \
                          IGRF_data.cxx \
                          JulianDate.cxx \
                          PointingHistory.cxx \
                          PointingInfo.cxx \
                          PointingTransform.cxx \
                          Quaternion.cxx \
                          SkyDir.cxx \
                          SkyProj.cxx \
                          SkyStat.cxx \
                          SolarSystem.cxx \
                          igrf_sub/igrf_sub.cxx \
                          jplephem/bary.cxx \
                          jplephem/dpleph.cxx

HD_CFLAGS               = ${HD_STD_CFLAGS}

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_CLEAN                = ${HD_LIBRARY_SRC_c:.c=.o} \
                          ${HD_LIBRARY_SRC_cxx:.cxx=.o}

include ${HD_STD_MAKEFILE}
