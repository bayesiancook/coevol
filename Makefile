all: _build/Makefile
	@cd _build ; make --no-print-directory

_build/Makefile:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build
	@rm -f cscope* *.dot tmp*
	@rm -rf data/tmp*
