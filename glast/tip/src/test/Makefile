HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_TEST_TASKS           = test_tip

HD_TEST_CXXTASK         = ${HD_TEST_TASKS}

HD_TEST_CXXTASK_SRC_cxx	= TestColumn.cxx \
                          TestExtensionData.cxx \
                          TestFileManager.cxx \
                          TestFileSummary.cxx \
                          TestFilter.cxx \
                          TestHarness.cxx \
                          TestImage.cxx \
                          TestInterpolation.cxx \
                          TestKeyRecord.cxx \
                          TestTable.cxx \
                          test_tip_main.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_TIP_LIBS} \
                          ${GT_FACILITIES_LIBS} \
                          ${ROOT_GUI_LIBS} \
                          ${HD_STD_LIBS} \
                          ${SYSLIBS}

include ${HD_STD_MAKEFILE}

# The following is to work around an unsolved problem in which the
# unit test seg faults when optimized.
test_tip_main.o: test_tip_main.cxx
	${CXX} -c -o $@ ${HD_CXXFLAGS} -g -O0 $?
