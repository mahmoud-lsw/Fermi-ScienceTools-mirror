HD_COMPONENT_NAME	= external

HD_COMPONENT_VERS	=

HD_SUBDIRS		= @h_components@ extFiles

HD_DISTCLEAN		= config.cache config.log config.status \
			  configure.lineno hd_install headas-setup hmakerc

HD_CLEAN		= ${HD_BLD_EXEC_PFX} ${HD_BLD_PFX} \
			  ../python/libpython*

HD_CTASK		= hd_install

HD_CTASK_SRC_c		= hd_install.c

HD_CFLAGS		= ${CFLAGS}

HD_INSTALL_TASKS	= hd_install

HD_INSTALL_SHELL_SCRIPTS	= hd_scanenv hmake

HD_INSTALL_EXTRA	= install-BUILD_DIR install-symlinks

include ${HD_STD_MAKEFILE}

distclean: hd-std-distclean
	cp -pf Makefile.in Makefile

# Install BUILD_DIR; tap dance needed because the built tools live in BUILD_DIR
#-------------------------------------------------------------------------------
install-BUILD_DIR:
	@${HD_INSTALL} ${HD_BUILD_DIR} ${HD_EXEC_PFX}/BUILD_DIR ${HD_CP_P}
#-------------------------------------------------------------------------------

# Specializations for lapack
#-------------------------------------------------------------------------------
subdir-lapack:
	cd ${LAPACK_DIR}; ${MAKE} ${HD_TARGET}

all-lapack default-lapack: build-lapack publish-lapack

build-lapack:
	${MAKE} subdir-lapack

publish-lapack:
	${MAKE} install-lapack HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-lapack:
	${MAKE} subdir-lapack HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB} includedir="${HD_INC}" \
		LIBDIR=${HD_LIB} INCDIR="${HD_INC}"

# Specializations for clhep
#-------------------------------------------------------------------------------
subdir-clhep:
	cd ${CLHEP_DIR}; ${MAKE} ${HD_TARGET}

all-clhep default-clhep: build-clhep publish-clhep

build-clhep:
	${MAKE} subdir-clhep

publish-clhep:
	${MAKE} install-clhep HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-clhep:
	${MAKE} subdir-clhep HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB}
#-------------------------------------------------------------------------------

# Specializations for cppunit
#-------------------------------------------------------------------------------
subdir-cppunit:
	cd ${CPPUNIT_DIR}; ${MAKE} ${HD_TARGET} ACLOCAL=: AMTAR=tar \
	  AUTOCONF=: AUTOHEADER=: AUTOMAKE=:

all-cppunit default-cppunit: build-cppunit publish-cppunit

build-cppunit:
	${MAKE} subdir-cppunit

publish-cppunit:
	${MAKE} install-cppunit HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-cppunit:
	${MAKE} subdir-cppunit HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB}

#------------------------------------------------------------------------------

# Specializations for f2c

#------------------------------------------------------------------------------

subdir-f2c:
	cd ${F2C_DIR}; ${MAKE} ${HD_TARGET}

all-f2c default-f2c: build-f2c publish-f2c

build-f2c:
	${MAKE} subdir-f2c

publish-f2c:
	${MAKE} install-f2c HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-f2c:
	${MAKE} subdir-f2c HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB} includedir="${HD_INC}/f2c" \
		LIBDIR=${HD_LIB} INCDIR="${HD_INC}/f2c"

#-------------------------------------------------------------------------------

# Specializations for fftw
#-------------------------------------------------------------------------------
subdir-fftw:
	cd ${FFTW_DIR}; ${MAKE} ${HD_TARGET}

all-fftw default-fftw: build-fftw publish-fftw

build-fftw:
	${MAKE} subdir-fftw

publish-fftw:
	${MAKE} install-fftw HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-fftw:
	${MAKE} subdir-fftw HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB} includedir="${HD_INC}/fftw"
#-------------------------------------------------------------------------------

# Specializations for python
#-------------------------------------------------------------------------------
subdir-python:
	cd ${PYTHON_DIR}; VPATH="" ${MAKE} ${HD_TARGET}

all-python default-python: build-python publish-python install-python

build-python:
	${MAKE} subdir-python

publish-python:
	${MAKE} install-python HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-python:
	${MAKE} subdir-python HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB}
#-------------------------------------------------------------------------------

# Specializations for root
#-------------------------------------------------------------------------------
subdir-root:
	@if [ "x${HD_TARGET}" = xall ]; then \
		cd ${ROOT_DIR}; ROOTSYS="${HD_EXEC_PFX}" ${MAKE}; \
	else \
		cd ${ROOT_DIR}; ROOTSYS="${HD_EXEC_PFX}" ${MAKE} ${HD_TARGET}; \
	fi

all-root default-root: build-root publish-root

build-root:
	${MAKE} subdir-root

publish-root:
	${MAKE} install-root HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-root:
	${MAKE} subdir-root HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB}
#-------------------------------------------------------------------------------

# Specializations for minuit2
#-------------------------------------------------------------------------------
subdir-minuit2:
	cd ${MINUIT2_DIR}; ${MAKE} ${HD_TARGET}

all-minuit2 default-minuit2: build-minuit2 publish-minuit2

build-minuit2:
	${MAKE} subdir-minuit2

publish-minuit2:
	${MAKE} install-minuit2 HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-minuit2:
	${MAKE} subdir-minuit2 HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB} includedir="${HD_INC}"
#-------------------------------------------------------------------------------

# Specializations for swig
#-------------------------------------------------------------------------------
subdir-swig:
	@if [ "x${HD_TARGET}" = xall ]; then \
		cd ${SWIG_DIR}; ${MAKE}; \
	else \
		cd ${SWIG_DIR}; ${MAKE} ${HD_TARGET}; \
	fi

all-swig default-swig: build-swig publish-swig

build-swig:
	${MAKE} subdir-swig

publish-swig:
	${MAKE} install-swig HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-swig:
	${MAKE} subdir-swig HD_TARGET="install" \
		prefix=${HD_EXEC_PFX} exec_prefix=${HD_EXEC_PFX} \
		libdir=${HD_LIB}
#-------------------------------------------------------------------------------

# Specializations for xerces
#-------------------------------------------------------------------------------
subdir-xerces:
	cd ${XERCES_DIR}; ${MAKE} ${HD_TARGET} XERCESCROOT="${HD_SRC}/xerces"

all-xerces default-xerces: build-xerces publish-xerces

build-xerces:
	${MAKE} subdir-xerces

publish-xerces:
	${MAKE} install-xerces HD_TOP_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_EXEC_PFX=${HD_BLD_EXEC_PFX} HD_CP_P="${HD_LN_S}"

install-xerces:
	${MAKE} subdir-xerces HD_TARGET="install" PREFIX=${HD_EXEC_PFX}
#-------------------------------------------------------------------------------
