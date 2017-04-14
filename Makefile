# ====================================
#                LISTS
# ====================================

SRC_FILES = $(shell find src -name "*.hpp") $(shell find src -name "*.cpp") $(shell find test -name "*.hpp") $(shell find test -name "*.cpp") $(shell find app -name "*.hpp") $(shell find app -name "*.cpp")
TMP_FILES = $(shell find . -name "tmp*")
.PHONY: cmake clean doc check format dot testdiffsel log perf report build-normal build-perf build-perffull build-coverage


# ====================================
#             COMPILATION
# ====================================
# Requires: cmake

all: cmake src/Eigen test/doctest.h
	@cd _build ; make --no-print-directory -j8

seq: cmake src/Eigen
	@cd _build ; make --no-print-directory

cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build doc/html
	@rm -f *.dot $(TMP_FILES)

fullclean: clean
	@rm -rf src/Eigen
	@rm -f test/doctest.h

src/Eigen:
	@wget --no-check-certificate http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
	@tar -xvf 3.3.3.tar.gz
	@cp -r eigen-eigen-67e894c6cd8f/Eigen src
	@rm -rf 3.3.3.tar.gz eigen-eigen-67e894c6cd8f

test/doctest.h:
	@wget --no-check-certificate https://raw.githubusercontent.com/onqtam/doctest/master/doctest/doctest.h
	@mv doctest.h test/


# ====================================
#              BUILD TYPES
# ====================================
# Requires: cmake

build-normal: clean cmake

build-perf:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=PERF ..

build-perffull:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=PERFFULL ..

build-coverage:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=COVERAGE ..

build-debug:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=DEBUG ..


# ====================================
#               TESTING
# ====================================
# Requires: graphviz (for dot)

test: all
	@cd _build ; make --no-print-directory test

testmove: all
	@_build/CustomDoubleMove

testdiffsel: all
	_build/diffsel data/c3c4/C4Amaranthaceaeshort.ali data/c3c4/C4Amaranthaceae.tree 3 5 tmp_diffsel_result clamp_MCMC 1 MS

testdiffsel2: all
	_build/diffsel data/ortho1to1_382g_26sp/FAM002112_1.fas.prank_codon.best.fas-gb.b.phylip data/ortho1to1_382g_26sp/FAM002112_1.fas.prank_codon.best.fas-gb.b.phylip_phyml_tree.txt 3 5 tmp_diffsel_result clamp_MCMC 1 MS

log:
	@less _build/Testing/Temporary/LastTest.log

dot: tmp.dot
	@dot -Tps $< -o tmp.ps
	@evince tmp.ps &


# ====================================
#        PERFORMANCE MEASUREMENT
# ====================================
# Requires: perf

perf: all
	@sudo bash -c 'echo "0" > /proc/sys/kernel/perf_event_paranoid' # nothing to see here :)
	@perf record _build/diffsel data/c3c4/C4Amaranthaceaeshort.ali data/c3c4/C4Amaranthaceae.tree 3 1 tmp_diffsel_result clamp_MCMC 1 MS

report: all
	@perf report | c++filt | less

report-save:
	@file=perf/report_`date +"%s"` ; echo "######################################\n# YOUR COMMENTS HERE:\n\n\n\n######################################\n" > "$$file" ; perf report | c++filt >> "$$file" ; nano "$$file"


# ====================================
#             CODE QUALITY
# ====================================
# Requires: clang-format, clang-check, clang-tidy

format:
	@clang-format-3.9 -i $(SRC_FILES)
# @clang-format -i $(SRC_FILES)

# WARNING: clang-tidy is not 100% reliable; use with caution!
# fix:
# 	@clang-tidy $(SRC_FILES) -checks=performance-* -fix -- -I src/ -std=gnu++11


# ====================================
#             DOCUMENTATION
# ====================================
# Requires: doxygen

doc:
	@doxygen Doxyfile
