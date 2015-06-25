#
#  Top Level Makefile for LAPACK
#  Version 3.4.1
#  April 2012
#

include make.inc

UNAME := $(shell uname)

ifneq ($(UNAME),Darwin)
all: lapack_install lib 
#blas_testing lapack_testing 

lib: blaslib lapacklib
	export LD_LIBRARY_PATH="${HD_BLD_EXEC_PFX}/lib:${LD_LIBRARY_PATH}"

distclean: cleanall

clean: cleanlib cleanblas_testing cleanall

lapack_install:
	( cd INSTALL; $(MAKE); ./testlsame; ./testslamch; ./testdlamch; \
	./testsecond; ./testdsecnd; ./testieee; ./testversion )

install: lib
	cp libblas${SHLIB_SUFFIX} $(LIBDIR)
	-ranlib $(LIBDIR)/libblas${SHLIB_SUFFIX}
	cp liblapack${SHLIB_SUFFIX} $(LIBDIR)
	-ranlib $(LIBDIR)/liblapack${SHLIB_SUFFIX}
#	cp libtmglib${SHLIB_SUFFIX} $(LIBDIR)
#	-ranlib $(LIBDIR)/libtmglib${SHLIB_SUFFIX}
#	-mkdir $(INCDIR)

blaslib:
	( cd BLAS/SRC; $(MAKE) )

lapacklib:	lapack_install
	( cd SRC; $(MAKE) )

else

all:

lib: blaslib lapacklib
		export LD_LIBRARY_PATH="${HD_BLD_EXEC_PFX}/lib:${LD_LIBRARY_PATH}"

distclean: cleanall

clean: cleanlib cleanblas_testing cleanall

endif

lapackelib: lapacklib
	( cd lapacke; $(MAKE) )

lapacke_example: lapackelib
	( cd lapacke/example; $(MAKE) )

variants:
	( cd SRC/VARIANTS ; $(MAKE))
	
tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

lapack_testing:	lib
	( cd TESTING ; $(MAKE) )
	./lapack_testing.py

variants_testing: lib variants
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/cholrl.a' ; \
	mv stest.out stest_cholrl.out ; mv dtest.out dtest_cholrl.out ; mv ctest.out ctest_cholrl.out ; mv ztest.out ztest_cholrl.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/choltop.a' ; \
	mv stest.out stest_choltop.out ; mv dtest.out dtest_choltop.out ; mv ctest.out ctest_choltop.out ; mv ztest.out ztest_choltop.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/lucr.a' ; \
	mv stest.out stest_lucr.out ; mv dtest.out dtest_lucr.out ; mv ctest.out ctest_lucr.out ; mv ztest.out ztest_lucr.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/lull.a' ; \
	mv stest.out stest_lull.out ; mv dtest.out dtest_lull.out ; mv ctest.out ctest_lull.out ; mv ztest.out ztest_lull.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/lurec.a' ; \
	mv stest.out stest_lurec.out ; mv dtest.out dtest_lurec.out ; mv ctest.out ctest_lurec.out ; mv ztest.out ztest_lurec.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/qrll.a' ; \
	mv stest.out stest_qrll.out ; mv dtest.out dtest_qrll.out ; mv ctest.out ctest_qrll.out ; mv ztest.out ztest_qrll.out )
		
blas_testing:
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 )
	( cd BLAS; ./xblat1s > sblat1.out    ; \
	           ./xblat1d > dblat1.out    ; \
	           ./xblat1c > cblat1.out    ; \
	           ./xblat1z > zblat1.out    ) 
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 )
	( cd BLAS; ./xblat2s < sblat2.in     ; \
	           ./xblat2d < dblat2.in     ; \
	           ./xblat2c < cblat2.in     ; \
	           ./xblat2z < zblat2.in     )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 )
	( cd BLAS; ./xblat3s < sblat3.in     ; \
	           ./xblat3d < dblat3.in     ; \
	           ./xblat3c < cblat3.in     ; \
	           ./xblat3z < zblat3.in     ) 
html:
	@echo "LAPACK HTML PAGES GENRATION with Doxygen"
	doxygen DOCS/Doxyfile
	@echo "=================="	
	@echo "LAPACK HTML PAGES GENRATED in DOCS/explore-html"
	@echo "Usage: open DOCS/explore-html/index.html"
	@echo "Online version available at http://www.netlib.org/lapack/explore-html/"
	@echo "=================="

man:
	@echo "LAPACK MAN PAGES GENRATION with Doxygen"
	doxygen DOCS/Doxyfile_man
	@echo "=================="
	@echo "LAPACK MAN PAGES GENRATED in DOCS/MAN"
	@echo "Set your MANPATH env variable accordingly"
	@echo "Usage: man dgetrf.f"
	@echo "=================="

cleanlib:
	( cd INSTALL; $(MAKE) clean )
	( cd BLAS/SRC; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
	( cd SRC/VARIANTS; $(MAKE) clean )
	

cleanblas_testing:	
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 clean )
	( cd BLAS; rm -f xblat* )


cleanall: cleanlib cleanblas_testing  
	rm -f *.a *.so *.dylib *${SHLIB_SUFFIX} TESTING/*.out INSTALL/test*  BLAS/*.out

