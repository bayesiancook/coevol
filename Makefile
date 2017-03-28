# ====================================
#                LISTS
# ====================================

SRC_FILES = $(shell find src -name "*.hpp") $(shell find src -name "*.cpp")
TMP_FILES = $(shell find . -name "tmp*")
.PHONY: cmake clean doc fix check format dot testmove testmove2 testgamma testdiffsel


# ====================================
#             COMPILATION
# ====================================
# Requires: cmake

all: cmake src/Eigen
	@cd _build ; make --no-print-directory -j8

cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build doc/html
	@rm -f *.dot $(TMP_FILES)

src/Eigen:
	@wget http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
	@tar -xvf 3.3.3.tar.gz
	@cp -r eigen-eigen-67e894c6cd8f/Eigen src
	@rm -rf 3.3.3.tar.gz eigen-eigen-67e894c6cd8f


# ====================================
#               TESTING
# ====================================
# Requires: graphviz (for dot)

test: all
	@cd _build ; make --no-print-directory test

testmove: all
	@_build/CustomSimpleMove

testmove2: all
	@_build/CustomDoubleMove

testgamma: all
	@_build/PoissonGamma data/test.data _build/test.out
	@less _build/test.out.trace

testdiffsel: all
	@_build/diffsel data/c3c4/C4Amaranthaceaeshort.ali data/c3c4/C4Amaranthaceae.tree 3 1 tmp_diffsel_result clamp_MCMC 1 MS

dot: tmp.dot
	@dot -Tps $< -o tmp.ps
	@evince tmp.ps &


# ====================================
#             CODE QUALITY
# ====================================
# Requires: clang-format, clang-check, clang-tidy

format:
	@clang-format -i $(SRC_FILES)

check:
	@clang-check $(SRC_FILES) -- -I src/ -std=gnu++11

# WARNING: clang-tidy is not 100% reliable; use with caution!
fix:
	@clang-tidy $(SRC_FILES) -checks=performance-* -fix -- -I src/ -std=gnu++11


# ====================================
#             DOCUMENTATION
# ====================================
# Requires: doxygen

doc:
	@doxygen Doxyfile
