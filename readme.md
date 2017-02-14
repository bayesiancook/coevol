# README #

## Compilation ##

`make`
: Builds mini-coevol (requires cmake); builds in the `_build` folder.

`make clean`
: Cleans the coevol directory; removes the `_build` folder among other things.

`make cmake`
: Runs cmake to build the makefile and `_build` folder but does not actually compile.


## Testing ##

`make test`
: Runs a two-step test. First, `poisson_gamma` is run and its output is displayed using `less`. Then `diffsel` is run; as it is very long, it should be interrupted after checking that it seems to run fine.

`make dot`
: Runs `poisson_gamma` in order to get a graph representation of the graphical model (requires `graphviz` and `evince`).


## Formatting ##

This code is formatted according to a specific style described in the `.clang-format` file (mostly Google-style with bigger indentation).

`make format`
: Automatically formats all the code according to the style file (this requires clang-format).

`make check`
: Runs clang-check on all the code (this requires clang-check). This can give non-trivial warnings that `gcc` usually doesn't catch.
