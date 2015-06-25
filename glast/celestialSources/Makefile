HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_SUBDIRS              = celestialSources \
                          src \
                          DarkMatter \
                          EarthPhenom \
						  eblAtten \
                          genericSources \
                          SpectObj \
                          GRB \
                          GRBobs \
                          microQuasar \
                          Pulsar \
                          GRBtemplate

HD_TEST_SUBDIRS			= src/test \
                          eblAtten \
                          genericSources \
                          GRB \
                          GRBobs \
                          GRBtemplate \
                          microQuasar

CERN_ROOT_SUBDIRS		= EarthPhenom \
						  GRB \
                          GRBobs \
                          GRBtemplate \
                          Pulsar \
                          SpectObj

include ${HD_STD_MAKEFILE}

# Custom target to prevent recursion to Root-specific packages as apppropriate.
subdir:
	exclude=`echo ${HD_SUBDIR} | sed -e 's:/.*::'`; \
	if [ "x${BUILD_WITHOUT_ROOT}" != "xyes" -o `echo ${CERN_ROOT_SUBDIRS} | grep -c $$exclude` -eq 0 ]; then \
	  ${HD_MAKE} ${HD_TARG_PFX}subdir; \
	fi
