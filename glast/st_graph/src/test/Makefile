HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_TEST_TASKS           = test_st_graph

HD_TEST_CXXTASK         = ${HD_TEST_TASKS}

HD_TEST_CXXTASK_SRC_cxx	= test_st_graph.cxx

HD_CXXFLAGS             = -fpermissive ${PYTHON_INC} ${HD_STD_CXXFLAGS}

HD_CXXLIBS              = ${HD_LFLAGS} \
                          ${GT_ST_GRAPH_LIBS} \
                          ${GT_ST_STREAM_LIBS} \
                          ${ROOT_GUI_LIBS} \
                          ${HD_STD_LIBS} \
                          ${SYSLIBS}

include ${HD_STD_MAKEFILE}
