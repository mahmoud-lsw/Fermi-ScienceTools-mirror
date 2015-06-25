HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_CXXTASK              = gtsrcid

HD_CXXTASK_SRC_cxx      = Catalogue.cxx \
                          Catalogue_fits.cxx \
                          Catalogue_id.cxx \
			  Catalogue_nr.cxx \
                          GHealpix.cxx \
                          GSkyDir.cxx \
                          Log.cxx \
                          Parameters.cxx \
                          sourceIdentify.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_CATALOGACCESS_LIBS} \
                          ${GT_STD_CXXLIBS}

HD_INSTALL_TASKS        = ${HD_CXXTASK}

HD_INSTALL_HELP         = gtsrcid.txt 

include ${HD_STD_MAKEFILE}
