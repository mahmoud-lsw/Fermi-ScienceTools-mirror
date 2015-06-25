HD_COMPONENT_NAME	= external

HD_COMPONENT_VERS	= 

all: default publish

default: build-numpy publish-numpy

publish: publish-numpy

install: install-numpy

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

# The LIBRARY_PATH environment variable is needed by Python packages
# with compiled C/C++ components, to find the Python libraries at
# compile time.

# N.B. Variables expressed as ${} are evaluated when this Makefile is
# processed. Variables expressed as $${} are evaluated when the _shell
# command_ is evaluated.

build-numpy:
	@echo 'BUILDING NUMPY...'; \
	export BLAS=${HD_BLD_LIB}; \
	export LAPACK=${HD_BLD_LIB}; \
	export PATH=${HD_TOP_BIN}:${HD_BIN}:${HD_BLD_BIN}:$${PATH}; \
	export LD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${LD_LIBRARY_PATH}; \
	export DYLD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${DYLD_LIBRARY_PATH}; \
	export LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${LIBRARY_PATH}; \
	python setup.py build --fcompiler=gnu95

publish-numpy:
	@echo 'PUBLISHING NUMPY...'; \
	export PATH=${HD_TOP_BIN}:${HD_BIN}:${HD_BLD_BIN}:$${PATH}; \
	export LD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${LD_LIBRARY_PATH}; \
	export DYLD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${DYLD_LIBRARY_PATH}; \
	python setup.py install

install-numpy:
	@echo 'INSTALLING NUMPY...'; \
	export PATH=${HD_TOP_BIN}:${HD_BIN}:${HD_BLD_BIN}:$${PATH}; \
	export LD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${LD_LIBRARY_PATH}; \
	export DYLD_LIBRARY_PATH=${HD_TOP_LIB}:${HD_LIB}:${HD_BLD_LIB}:$${DYLD_LIBRARY_PATH}; \
	python setup.py install
