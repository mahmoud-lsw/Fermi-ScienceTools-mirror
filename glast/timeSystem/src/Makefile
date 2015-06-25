HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = timeSystem

HD_LIBRARY_SRC_c        = axBary.c \
                          bary.c \
                          clock.c \
                          ctatv.c \
                          dpleph.c \
                          glastscorbit.c \
                          phaseHist.c \
                          scorbit.c \
                          xCC.c

HD_LIBRARY_SRC_cxx      = AbsoluteTime.cxx \
                          BaryTimeComputer.cxx \
                          CalendarFormat.cxx \
                          Duration.cxx \
                          ElapsedTime.cxx \
                          EventTimeHandler.cxx \
                          GlastTimeHandler.cxx \
                          IntFracUtility.cxx \
                          MjdFormat.cxx \
                          PulsarTestApp.cxx \
                          SourcePosition.cxx \
                          TimeCorrectorApp.cxx \
                          TimeInterval.cxx \
                          TimeSystem.cxx

HD_SUBDIRS              = gtbary

HD_CFLAGS               = ${HD_STD_CFLAGS}

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

include ${HD_STD_MAKEFILE}
