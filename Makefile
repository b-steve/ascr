ALL:
	make prepare
	make rcpp
	make roxygen
	make check
	make install

prepare:
	rm -rfv man
	rm -fv NAMESPACE src/*.o src/RcppExports.cpp src/admbsecr.so src/symbols.rds R/RcppExports.R

rcpp:
	R --slave -e "library(Rcpp); compileAttributes()"

roxygen:
	R --slave -e "library(roxygen2); roxygenise('~/admbsecr/')"

check:
	R CMD check .

install:
	R CMD INSTALL .

clean:
	rm -rfv ..Rcheck/
	rm -rfv src/*.o src/*.so src/*.rds


