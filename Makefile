all: _build/Makefile
	@cd _build ; make --no-print-directory

_build/Makefile:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build
	@rm -f cscope* *.dot tmp* src/cscope*
	@rm -rf data/tmp*

test: all
	@_build/PoissonGamma data/test.data _build/test.out
	@less _build/test.out.trace
	@_build/DiffSel data/c3c4/C4Amaranthaceaeshort.ali data/c3c4/C4Amaranthaceae.tree 3 1 tmp_diffsel_result clamp_MCMC 1 MS
