# README #

## Getting started ##

### Compilation instructions ###
To build __coevol__/__diffsel__, simply run:

`make`

from the root of the __coevol__ folder (ie, the folder that contains this readme).
If you get errors, check the requirements section below.

To run a series of simple tests to check that __coevol__ is compiled correctly, run `make test`.
You can also run `make testdiffsel` to run the diffsel application on the c3c4 data in order to
check that everything seems to run fine (you'll probably want to interrupt it after a few iterations).


### Requirements ###
This version of __coevol__ requires a decently recent version of __cmake__ and an internet connection (the ability to run `wget` commands).

If you get a message saying something resembling `cmake command not found`, please install __cmake__ (`sudo apt-get install cmake` on debian/ubuntu).

If you get a message that looks like this:

```
CMake Error at CMakeLists.txt:1 (cmake_minimum_required):
    CMake 3.1.0 or higher is required.  You are running version 2.5.2
```
    
then it is possible that the default __cmake__ on your machine is too old. Please send an email to `vincent.lanore (at) univ-lyon1.fr`
and include the output of `cmake --version` and the output of `make clean; make`.


### Running diffsel ###
After compilation, there should be a `diffsel` executable in the `_build/` folder (itself located at the root folder of __coevol__).
This __diffsel__ executable expects the following arguments (in that order):
* an alignment file (supports the phylip format);
* a tree file (supports the newick format); branch labels should be integers (between 0 and P-1), specifying a partition of the set of branches into P subsets
* the number of conditions K (should be at most equal to number of partitions P defined by the tree; can be less than P, in which case all partitions with index greater than K will be allocated to last condition)
* how many iterations to perform before writing to disk (eg, value 5 will save to disk every 5 iterations);
* the name of the run (used to name the output files);
* either `clamp_MCMC`, `clamp_MCMC_var`, or `unclamp`:
	* `clamp_MCMC`: uniform distribution of amino-acid profiles across sites, variance of differential selection effects  estimated from the data (recommended by default) 
	* `clamp_MCMC_var`: uniform distribution of amino-acid profiles across sites, variance of differential selection effects fixed to `1` (faster, but less easily justifiable);	
	* `unclamp`: all hyperparameters estimated from data (currently not recommended)
* value of the conjugate parameter 
	* `1` : conjugate sampling activated (recommended by default)
	* `0` : conjugate sampling inactivated;
* either `SR` (square root) of `MS` (mutation selection) model for the relation between fitness parameters and substitution rates between codons (`MS` recommended by default).

For example, this is the command run by `make testdiffsel`:

```bash
_build/diffsel data/c3c4/C4Amaranthaceaeshort.ali data/c3c4/C4Amaranthaceae.tree 3 5 tmp_diffsel_result clamp_MCMC 1 MS
```


## Makefile commands ##

### Compilation ###
`make`
: Builds __coevol__ in parallel (requires `cmake`); builds in the `_build` folder.

`make seq`
: Builds __coevol__ sequentially. Can be useful to catch errors one by one or to avoid occupying too many resources.

`make clean`
: Cleans the coevol directory; removes the `_build` folder among other things.

`make cmake`
: Runs cmake to build the makefile and `_build` folder but does not actually compile.


### Testing ###
`make test`
: Runs all tests registered to `cmake`, ie, all tests from the `test/` folder.

`make dot`
: Displays a graph representation of the graphical model of the last application that invoked `getDot` (requires `graphviz` and `evince`). You can run `_build/PoissonGamma` and then`make dot` to get an idea of what this does.


### Code quality ###
This code is formatted according to a specific style described in the `.clang-format` file (mostly Google-style with bigger indentation).

`make format`
: Automatically formats all the code according to the style file (this requires `clang-format`).

`make check`
: Runs clang-check on all the code (this requires `clang-check`). This can give non-trivial warnings that `gcc` usually doesn't catch.


### Documentation ###
__coevol__ provides a basic Doxygen documentation.

`make doc`
: Generates the Doxygen documentation (requires `doxygen`) in `doc/html`.
