ALL:
	make compile
	make prepare
	make rcpp
	make roxygen
	make check
	make install

compile: inst/ADMB/src/densfuns.cpp inst/ADMB/src/detfuns.cpp inst/ADMB/src/secr.tpl
	cd inst/ADMB/src; admb -O secr.tpl
	cd inst/ADMB/src; rm -rfv secr.cpp secr.htp secr.o secr.obj
	cd inst/ADMB/bin/linux; rm -rfv secr_diff; mv ../../src/secr ./secr_diff

prepare:
	rm -rfv man
	rm -fv NAMESPACE src/*.o src/RcppExports.cpp src/admbsecr.so src/symbols.rds R/RcppExports.R

rcpp:
	R --slave -e "library(Rcpp); compileAttributes()"

roxygen:
	R --slave -e "library(roxygen2); roxygenise('/scratch/admbsecr/')"

check:
	R CMD check .

install:
	R CMD INSTALL .

clean:
	rm -rfv ..Rcheck/
	rm -rfv src/*.o src/*.so src/*.rds
	rm -rfv src-i386/ src-x64



