# ====================================
#              VARIABLES
# ====================================

SRC_FILES = $(shell ls include/**/*.hpp) $(shell ls src/**/*.cpp) $(shell ls src/*.cpp)


# ====================================
#             COMPILATION
# ====================================
# Requires: cmake

all: cmake
	@cd _build ; make --no-print-directory

cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build
	@rm -f *.dot tmp*
	@rm -rf data/tmp*


# ====================================
#               TESTING
# ====================================
# Requires: graphviz (for dot)

test: all
	@_build/poisson_gamma data/test.data _build/test.out
	@less _build/test.out.trace
	@_build/diffsel data/c3c4/C4Amaranthaceaeshort.ali data/c3c4/C4Amaranthaceae.tree 3 1 tmp_diffsel_result clamp_MCMC 1 MS

dot: all
	@_build/poisson_gamma data/test.data _build/test.out
	@dot -Tps tmp.dot -o tmp.ps
	@evince tmp.ps &


# ====================================
#             CODE QUALITY
# ====================================
# Requires: clang-format, clang-check, clang-tidy

format:
	@clang-format -i $(SRC_FILES)

check:
	@clang-check $(SRC_FILES) -- -I include/ -std=gnu++11

# WARNING: clang-tidy is not 100% reliable; use with caution!
fix:
	@clang-tidy $(SRC_FILES) -checks=performance-* -fix -- -I include/ -std=gnu++11
