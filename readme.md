# README #

## Compilation ##

`make`
: Builds __`mini-coevol`__ in parallel (requires `cmake`); builds in the `_build` folder.

`make seq`
: Builds __`mini-coevol`__ sequentially. Can be useful to catch errors one by one or to avoid occupying too many resources.

`make clean`
: Cleans the coevol directory; removes the `_build` folder among other things.

`make cmake`
: Runs cmake to build the makefile and `_build` folder but does not actually compile.


## Testing ##

`make test`
: Runs all tests registered to `cmake`, ie, all tests from the `test/` folder.

`make dot`
: Displays a graph representation of the graphical model of the last application that invoked `getDot` (requires `graphviz` and `evince`). You can run `make testmove` and then`make dot` to get an idea of what this does.


## Code quality ##

This code is formatted according to a specific style described in the `.clang-format` file (mostly Google-style with bigger indentation).

`make format`
: Automatically formats all the code according to the style file (this requires `clang-format`).

`make check`
: Runs clang-check on all the code (this requires `clang-check`). This can give non-trivial warnings that `gcc` usually doesn't catch.


## Documentation ##

__`mini-coevol`__ provides a basic Doxygen documentation.

`make doc`
: Generates the Doxygen documentation (requires `doxygen`) in `doc/html`.
