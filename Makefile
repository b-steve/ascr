ALL:
	make compile
	make prepare
	make rcpp
	make roxygen
	make clean
	make build
	make check
	make install

docs:
	make compile
	make prepare
	make rcpp
	make roxygen

compile: inst/ADMB/src/densfuns.cpp inst/ADMB/src/detfuns.cpp inst/ADMB/src/secr.tpl
	if [ $(shell hostname) == "heton" ]; then cd inst/ADMB/src; admb -O secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/linux; mv ../../src/secr ./secr; fi
	if [ $(shell hostname) == "albatross.mcs.st-and.ac.uk" ]; then cd inst/ADMB/src; admb -O secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/mac; rm -rfv secr; mv ../../src/secr ./secr; fi
	if [ $(shell hostname) == "morten-win" ]; then export OLDPATH=$$PATH; export PATH=/c/"Program Files (x86)"/ADMB/bin:/c/"Program Files (x86)"/ADMB/utilities:/c/"Program Files (x86)"/ADMB/utilities/mingw64/bin:/c/admb/bin:$$PATH; echo $$PATH; cd inst/ADMB/src; cmd //c admb -f secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/windows; rm -rfv secr.exe; mv ../../src/secr.exe ./secr.exe; export PATH=$$OLDPATH; fi
	if [ $(shell hostname) == "midge" ]; then export OLDPATH=$$PATH; echo $$OLDPATH; export PATH=/c/admb/bin:/c/MinGW/bin:$$PATH; cd inst/ADMB/src; admb -O secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/windows; rm -rfv secr.exe; mv ../../src/secr.exe ./secr.exe; export PATH=$$OLDPATH; fi

prepare:
	rm -rfv man
	rm -fv NAMESPACE src/*.o src/RcppExports.cpp src/admbsecr.so src/symbols.rds R/RcppExports.R

rcpp:
	R --slave -e "library(Rcpp); compileAttributes()"

roxygen:
	R --slave -e "library(roxygen2); roxygenise('.')"

build:
	if [ $(shell hostname) == "morten-win" ]; then rm -rfv .Rbuildignore; fi
	if [ $(shell hostname) == "midge" ]; then rm -rfv .Rbuildignore; fi
	R CMD build --resave-data .
	if [ $(shell hostname) == "morten-win" ]; then mkdir .Rbuildignore; fi
	if [ $(shell hostname) == "midge" ]; then mkdir .Rbuildignore; fi
	mv admbsecr_1.1.1.tar.gz .Rbuildignore/

check:
	R CMD check .Rbuildignore/admbsecr_1.1.1.tar.gz --no-tests

install:
	R CMD INSTALL .Rbuildignore/admbsecr_1.1.1.tar.gz --install-tests

pdf:
	R CMD Rd2pdf --pdf . &
	rm ..pdf

clean:
	rm -rfv ..Rcheck/ ..pdf
	rm -rfv src/*.o src/*.so src/*.rds
	rm -rfv src-i386/ src-x64/
	rm -rfv admbsecr.Rcheck

