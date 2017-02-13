SRC_FILES = $(shell ls include/**/*.hpp) $(shell ls src/**/*.cpp) $(shell ls src/*.cpp)

all: cmake
	@cd _build ; make --no-print-directory

cmake : _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build
	@rm -f *.dot tmp*
	@rm -rf data/tmp*

test: all
	@_build/poisson_gamma data/test.data _build/test.out
	@less _build/test.out.trace
	@_build/diffsel data/c3c4/C4Amaranthaceaeshort.ali data/c3c4/C4Amaranthaceae.tree 3 1 tmp_diffsel_result clamp_MCMC 1 MS

dot: all
	@_build/poisson_gamma data/test.data _build/test.out
	@dot -Tps tmp.dot -o tmp.ps
	@evince tmp.ps &

format:
	@clang-format -i $(SRC_FILES)
