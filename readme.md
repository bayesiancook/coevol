# README #

## Getting started ##

### Compilation instructions ###

To build __coevol__/__diffsel__, simply run:

`make`

from the root of the __coevol__ folder.
If you get errors, check the requirements section below.


### Requirements ###

This version of __coevol__ requires a decently recent version of __cmake__ and an internet connection (the ability to run `wget` commands).

If you get a message saying something like `cmake command not found`, please install __cmake__ (`sudo apt-get install cmake` on debian/ubuntu).

If you get a message that looks like this:

```
CMake Error at CMakeLists.txt:1 (cmake_minimum_required):
    CMake 3.1.0 or higher is required.  You are running version 2.5.2```
    
then it is possible that the default __cmake__ on your machine is too old. Please send an email to `vincent.lanore@univ-lyon1.fr`
and include the output of `cmake --version` and the output of `make clean; make`.


## Makefile commands ##

### Compilation ###

`make`
: Builds __`mini-coevol`__ in parallel (requires `cmake`); builds in the `_build` folder.

`make seq`
: Builds __`mini-coevol`__ sequentially. Can be useful to catch errors one by one or to avoid occupying too many resources.

`make clean`
: Cleans the coevol directory; removes the `_build` folder among other things.

`make cmake`
: Runs cmake to build the makefile and `_build` folder but does not actually compile.


### Testing ###

`make test`
: Runs all tests registered to `cmake`, ie, all tests from the `test/` folder.

`make dot`
: Displays a graph representation of the graphical model of the last application that invoked `getDot` (requires `graphviz` and `evince`). You can run `make testmove` and then`make dot` to get an idea of what this does.


### Code quality ###

This code is formatted according to a specific style described in the `.clang-format` file (mostly Google-style with bigger indentation).

`make format`
: Automatically formats all the code according to the style file (this requires `clang-format`).

`make check`
: Runs clang-check on all the code (this requires `clang-check`). This can give non-trivial warnings that `gcc` usually doesn't catch.


### Documentation ###

__`mini-coevol`__ provides a basic Doxygen documentation.

`make doc`
: Generates the Doxygen documentation (requires `doxygen`) in `doc/html`.
