HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_XML_LIBS} \
                          ${GT_FACILITIES_LIBS}

HD_TEST_TASKS           = entity_test \
                          test_xmlBase \
                          test_IFile \
                          test_mem \
                          test_parse \
                          test_write

HD_CLEAN                = entity_test${BSUF} \
                          test_xmlBase${BSUF} \
                          test_IFile${BSUF} \
                          test_mem${BSUF} \
                          test_write${BSUF}

include ${HD_STD_MAKEFILE}

test:
	${HD_MAKE} entity_test HD_TEST_CXXTASK=entity_test HD_TEST_CXXTASK_SRC_cxx=entity_test.cxx
	${HD_MAKE} test_xmlBase HD_TEST_CXXTASK=test_xmlBase HD_TEST_CXXTASK_SRC_cxx=main.cxx
	${HD_MAKE} test_IFile HD_TEST_CXXTASK=test_IFile HD_TEST_CXXTASK_SRC_cxx=test_IFile.cxx
	${HD_MAKE} test_mem HD_TEST_CXXTASK=test_mem HD_TEST_CXXTASK_SRC_cxx=test_mem.cxx
	${HD_MAKE} test_write HD_TEST_CXXTASK=test_write HD_TEST_CXXTASK_SRC_cxx=test_write.cxx
