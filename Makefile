SHELL := /bin/bash
ALL:
	make compile
	make nocompile

nocompile:
	make prepare
	make roxygen
	make rcpp
	make roxygen
	make clean
	make build
	make check
	make install

docs:
	make compile
	make prepare
	make roxygen
	make rcpp

compile: inst/ADMB/src/densfuns.cpp inst/ADMB/src/detfuns.cpp inst/ADMB/src/secr.tpl
	if [ $(shell hostname) == "sc-stat-369411" ]; then cd inst/ADMB/src; LDFLAGS=-static admb -f secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/linux; rm -rfv secr; mv ../../src/secr ./secr; fi
	if [ $(shell hostname) == "darwen" ]; then cd inst/ADMB/src; LDFLAGS=-static admb -f secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/linux; rm -rfv secr; mv ../../src/secr ./secr; fi
	if [ $(shell hostname) == "fygwyd" ]; then cd inst/ADMB/src; LDFLAGS=-static admb -f secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/linux; rm -rfv secr; mv ../../src/secr ./secr; fi
	if [ $(shell hostname) == "darwen-win" ]; then export OLDPATH=$$PATH; export PATH=/c/admb:/c/admb/utilities/mingw64/bin:$$PATH; echo $$PATH; cd inst/ADMB/src; cmd //c admb -f secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/windows; rm -rfv secr.exe; mv ../../src/secr.exe ./secr.exe; export PATH=$$OLDPATH; fi
	if [ $(shell hostname) == "fygwyd-win" ]; then cd inst/ADMB/src; cmd //c admb -f secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/windows; rm -rfv secr.exe; mv ../../src/secr.exe ./secr.exe; fi
	if [ $(shell hostname) == "macbook-pro-di-franchini-filippo-2.home" ]; then export OLDPATH=$$PATH; export PATH=/Documents/admb/admb/build/dist/bin/:/Documents/admb/admb/build/dist/utilities/:$$PATH; echo $$PATH; cd inst/ADMB/src; ~/Documents/admb/admb/admb -f secr.tpl; rm -rfv secr.cpp secr.htp secr.o secr.obj; cd ../bin/mac; rm -rfv secr; mv ../../src/secr ./secr; fi

prepare:
	rm -rfv man
	rm -fv NAMESPACE src/*.o src/RcppExports.cpp src/ascr.so src/symbols.rds R/RcppExports.R

rcpp:
	export str="library(Rcpp);compileAttributes()"; R --slave -e $$str

roxygen:
	export str="library(roxygen2);roxygenise('.')"; R --slave -e $$str

build:
	R CMD build --resave-data .
	mkdir -p package-build
	mv ascr_2.2.0.tar.gz package-build/

check:
	R CMD check package-build/ascr_2.2.0.tar.gz --no-tests --no-examples

install:
	R CMD INSTALL package-build/ascr_2.2.0.tar.gz --install-tests

pdf:
	R CMD Rd2pdf --pdf . &
	rm ..pdf

clean:
	rm -rfv ..Rcheck/ ..pdf
	rm -rfv src/*.o src/*.so src/*.rds
	rm -rfv src-i386/ src-x64/
	rm -rfv ascr.Rcheck

