HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = xmlBase

HD_LIBRARY_SRC_cxx      = Dom.cxx \
                          EResolver.cxx \
                          IFile.cxx \
                          XmlErrorHandler.cxx \
                          XmlParser.cxx \
                          docMan/DocMan.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

HD_CLEAN                = ${HD_LIBRARY_SRC_cxx:.cxx=.o}

include ${HD_STD_MAKEFILE}
