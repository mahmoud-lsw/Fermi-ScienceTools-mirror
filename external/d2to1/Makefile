HD_COMPONENT_NAME	= external

HD_COMPONENT_VERS	= 

all: publish

#default: build-d2to1

publish: publish-d2to1

install: install-d2to1

include ${HD_STD_MAKEFILE}

# The special targets below are needed to properly wrap the 'python
# setup.py' process in the hmake structure. The 'python setup.py'
# standard is specific to the Python world, and does not play well
# with other build systems, IMHO.

# The code below which sets PATH ensures that a HEADAS version of
# Python is found even if only a single package is being built,
# e.g. rebuilding a single package after complete installation. The
# same reasoning applies for the library paths, since the python
# binary must be able to find the python shared library.

# N.B. Variables expressed as ${} are evaluated when this Makefile is
# processed. Variables expressed as $${} are evaluated when the _shell
# command_ is evaluated.

build-d2to1:
	@echo 'BUILDING D2TO1...'; \
	export PATH=${HD_TOP_BIN}:${HD_BIN}:${HD_BLD_BIN}:$${PATH}; \
	export LD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${LD_LIBRARY_PATH}; \
	export DYLD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${DYLD_LIBRARY_PATH}; \
	python setup.py build

publish-d2to1:
	@echo 'PUBLISHING D2TO1...'; \
	export PATH=${HD_TOP_BIN}:${HD_BIN}:${HD_BLD_BIN}:$${PATH}; \
	export LD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${LD_LIBRARY_PATH}; \
	export DYLD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${DYLD_LIBRARY_PATH}; \
	python setup.py install

install-d2to1:
	@echo 'INSTALLING D2TO1...'; \
	export PATH=${HD_TOP_BIN}:${HD_BIN}:${HD_BLD_BIN}:$${PATH}; \
	export LD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${LD_LIBRARY_PATH}; \
	export DYLD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${DYLD_LIBRARY_PATH}; \
	python setup.py install
