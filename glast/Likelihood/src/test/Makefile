HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_TEST_TASKS           = test_Likelihood

HD_TEST_CXXTASK         = ${HD_TEST_TASKS}

HD_TEST_CXXTASK_SRC_cxx	= test.cxx \
                          SourceData.cxx \
                          XmlDiff.cxx

#------------------------------------------------------------------------------

# <HACK>

# This hack is only needed to get test.cxx to build under hmake. The
# default compiler switche variable HD_STD_CXXFLAGS contains
# '--pedantic', which causes code in test.cxx to not compile. So I
# strip out that switch for now. Eventually, the code in test.cxx
# should be fixed to address this issue.

# N.B. This line is the ONLY line needed to implement this hack. When
# test.cxx iis fixed, this line may be removed.

# Eric Winter, 21 June 2010

HD_STD_CXXFLAGS         = `echo $$CXXOPT $$CXXFLAGS $$HD_IFLAGS $$DEFS | sed -e "s/--pedantic//"`

# </HACK>

#------------------------------------------------------------------------------

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS} -fpermissive

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_LIKELIHOOD_LIBS} \
			  -leblAtten \
                          ${GT_MAP_TOOLS_LIBS} \
                          ${GT_DATASUBSELECTOR_LIBS} \
                          ${GT_EVTBIN_LIBS} \
                          ${GT_OPT_LIBS} \
                          ${GT_HEALPIX_LIBS} \
                          -lirfLoader \
                          -lirfInterface \
                          -lirfUtil \
                          -llatResponse \
                          -ldc1aResponse \
                          -ltestResponse \
                          ${GT_XML_LIBS} \
                          ${GT_FFT_LIBS} \
                          ${GT_CPPUNIT_LIBS} \
                          ${GT_STD_CXXLIBS}

include ${HD_STD_MAKEFILE}

#Do not change this.  It is needed to work around the behavior
#of make when using the test target and a file named test.cxx
test:
	${HD_MAKE} ${HD_TEST_CXXTASK}
