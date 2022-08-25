dune-nonlinopt: General-Purpose Unconstrained Nonlinear Optimization
====================================================================

The package dune-nonlinopt is a flexible, efficient and robust
toolkit for unconstrained nonlinear optimization. In particular, it
provides implementations of the following optimization schemes:

1. The limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)
    method [(ref)][1]
2. The nonlinear generalized minimal residual (N-GMRES) method
    [(ref)][2]
3. Nonlinear conjugate gradients (NCG) in five different variants,
    namely Hestenes-Stiefel, Fletcher-Reeves, Polak-Ribière-Polyak,
    Dai-Yuan and Hager-Zhang [(ref)][3]
4. The Barzilai-Borwein method [(ref)][4] and good old steepest descent
5. Flexible combinations of these methods, e.g., L-BFGS in conjugate
    gradients mode (SCG, [ref][5]), or N-GMRES with conjugate gradients
    as preconditioner

The default line search employed by dune-nonlinopt is a
reimplementation of Hager and Zhang's `CG_DESCENT` line search,
combined with nonmonotone ``relaxed'' Wolfe conditions based on the
average of previous function values [(ref)][6] and the approximate Wolfe
conditions [(ref)][7]. Alternative line searches (cubic, quadratic,
backtracking) and line search termination criteria (e.g., several
variants of the Wolfe conditions) are also available.

Here is a list of the main features offered by dune-nonlinopt,
many of which are not found in other large optimization packages:

* Support for arbitrary scalar and vector types, not just single
    / double precision, but also arbitrary precision arithmetic,
    or optimization on high-dimensional parameter fields with
    flexible lower-dimensional representations (e.g., local refinement).
* Full support for large-scale parallel computations, e.g., within
    the [DUNE framework][8]. The user only has to make sure that the
    arithmetic operations of the chosen vector class, like scalar
    products, are parallelized.
* A general-purpose preconditioner interface, which makes it
    possible to employ preconditioned variants of the aforementioned
    methods.
* A flexible plugin framework, with the possibility to exchange
    virtually all subalgorithms with user-defined versions, e.g.,
    custom line searches and termination criteria.

In contrast, many existing optimization codes use hard-coded
double precision arrays for the representation of parameter vectors
and explicit low-level operations on them, which makes it difficult
or unfeasable to exchange the underlying data type or perform
parallel computations. Few software packages offer direct support for
preconditioning, and those that offer configurable subalgorithms
are often restricted to a small pool of predefined variants.


Copying and Contributing
========================

dune-nonlinopt is free software under the BSD 3-clause license, and its
documentation is licensed under a Creative Commons license. See
[LICENSE.md](./LICENSE.md) for details.

If you would like to contribute, please feel free to open a merge
request. It would of course be nice if your code would follow the
style conventions already present in the source, but there are
currently no strict guidelines beyond that. You are invited to use
the bug tracker of this repository if you find any errors or mistakes,
and you may also contact the author directly
(ole.klein@iwr.uni-heidelberg.de) if there is something you would like
to discuss in detail.


Dependencies
============

The package does not have any hard dependencies, except for a
not-too-ancient C++ compiler (C++17 and above), a reasonably new
CMake installation (3.10 and above), and a version of make.
dune-nonlinopt is header-only, so technically one could even use it
with just the C++ compiler and no automatic build support.

Optional Dependencies
---------------------

* The N-GMRES method and those derived from it require the
    [Eigen3 package][9], since its QR decomposition is used in solving
    the local least squares problems.
* The optional support for configuration files and string-based
    configuration requires the [dune-common][10] package.
* The optional parser-based problem definition requires
    [muParser][11]. This can be used to specify comparatively simple
    optimization problems through mathematical expressions instead
    of having to implement a C++ class, see the examples section.
* Building the automatic documentation requires Doxygen.


Installation
============

dune-nonlinopt supports three different ways of installation:
1. Installation as a submodule of the [DUNE framework][8]
2. Using CMake for configuration and installation
3. A barebones header-only approach

Integration into DUNE (Default)
-------------------------------

dune-nonlinopt can be built and used as a DUNE module, with the DUNE
installer handling dependencies. The only other required module is
[dune-common][10]. Simply place the source next to your other DUNE
modules, and let the `dunecontrol` script take care of things. Here
is an example options file:

```bash
# subdirectory to use as build-directory
BUILDDIR="$HOME/dune/build"
 
# paths to external software in non-default locations
CMAKE_PREFIX_PATH="$HOME/software"
 
GXX_WARNING_OPTS="-Wall"
 
GXX_OPTS="-march=native -g -O2 -std=c++20"
 
# set flags for CMake:
# compiler flags, 3rd party software, grid extensions
CMAKE_FLAGS=" \
-DCMAKE_Fortran_COMPILER='gfortran-10.2.0' \
-DCMAKE_C_COMPILER='gcc-10.2.0' \
-DCMAKE_CXX_COMPILER='g++-10.2.0' \
-DCMAKE_CXX_FLAGS='$GXX_WARNING_OPTS $GXX_OPTS' \
-DDUNE_SYMLINK_TO_SOURCE_TREE=1 \
"
```

Most of these options are not really necessary, they are there to
show you how DUNE can be configured. Have a look at the DUNE
webpage if you want to know more, and make sure to either
adapt or remove the compiler specifications above. Here is how you
build the project using an options file called `myopts`:

```bash
dune-common/bin/dunecontrol --opts=myopts --module=dune-nonlinopt
```

Afterwards, you can run `make build_tests` in the build directory
of dune-nonlinopt to build the automatic tests, if desired, and
then execute them with `ctest`. These test problems are a subset
of the ones found in [(ref)][12]. The automatic documentation can be
built by issuing `make doc` if Doxygen was found. Open
`doc/html/index.html` in a browser of your choice to read it.

Building with CMake
-------------------

You can also build the package directly with CMake. Create a build
directory, either as a subfolder within the source directory tree,
or in some other location, and then call CMake to configure the
project and build it with make:

```bash
# create build directory
mkdir build
cd build

# run CMake for the build directory
# replace with path to source if using custom location
cmake ..

# optional: modify configuration (see below)
ccmake .

# build the project
make

# optional: install header files (see below)
make install
```

CMake will try to detect optional dependencies and generate an
initial configuration. Calling `ccmake .` within the build directory,
you can modify this configuration to your liking. Note that some
distros don't install ccmake with make, but provide it as a separate
package instead. You can then either install ccmake directly, or set
the following options noninteractively by using the `-D` commandline
option when calling cmake as described above.

Enable / disable optional parts of the project:
* `BUILD_DOCUMENTATION` toggles whether or not Doxygen documentation
    is generated within the `doc/` subfolder. If you build documentation,
    you may access it by opening `doc/html/index.html` in a browser of
    your choice.
* `BUILD_EXAMPLES` toggles whether the toy problems within the
    `examples` folder are built. These may be a good starting point for
    your own projects. Note that the problemparser example requires
    additional software to be installed, namely the muParser library.
* `BUILD_TESTING` enables the compilation of tests within the `test/`
    directory. You can use `ctest` to run these tests and sanity checks.
    Note that only a subset of the test suite will be built and run
    this way, since some of the tests require string-based configuration,
    and therefore dune-common.

Configure the CMake build process:
* `CMAKE_BUILD_TYPE` is used to choose compilation options. The default
    of CMake is an unoptimized `Debug` build, which is helpful when you
    encounter bugs while writing your own application. Be sure to switch
    this over to `Release` to turn optimization on once everything works
    as expected.
* `CMAKE_INSTALL_PREFIX` specifies the directory where the project should
    be installed if you type `make install`.
* `CMAKE_PREFIX_PATH` can be used to guide the build system when you
    have installed software in non-standard locations and would like to
    make it available.

Using `make install`, the header files in the `dune/nonlinopt/` directory
can be copied to the subdirectory `include/dune-nonlinopt/` of the path
specified by the `CMAKE_INSTALL_PREFIX` variable. This step is entirely
optional. Note that the example and test executables are not installed by
this command.


Header Only Approach
--------------------

Except for the examples and automatic tests contained in their
respective directories, dune-nonlinopt does not contain any libraries
or executables. It is therefore straight-forward to ``install''
the package by copying the contents of the `src` folder to a
destination that is searched by the compiler for header files,
or to a folder that is then made known to the compiler through
an appropriate flag, e.g., `-I/path/to/headers` or similar.

As a special case, the repository can be incorporated into
another project as a git submodule.


Examples
========

The `examples` folder contains toy problems that illustrate how 
dune-nonlinopt can be used to solve optimization problems.

Included are:

* `problemparser`: Parses a user-supplied mathematical expression
    defining a cost function, and solves the resulting minimization
    problem for a given starting position. Requires muParser.
* `rosenbrock`: Demonstrates different ways of creating and
    configuring solver objects. Different driver functions show how
    parts of the algorithms can be configured, including parameter
    values and potentially user-supplied custom classes.
    String-based configuration is also shown, but building and
    running this part requires dune-common.


Quick Start Instructions
========================

This is a description of the general workflow when solving an
optimization problem with dune-nonlinopt. Most of these steps
are also demonstrated within the example files mentioned above.
Have a look at the automatic documentation if you would like
to know more about the class information and definitions.

Specifying the Optimization Problem
-----------------------------------

* If the objective function is a simple mathematical expression,
    consider using the `problemparser` example program, or basing
    your code on the `ParsedProblem` C++ class. This way, you don't
    need to write any custom methods.
* In other cases, derive from `ProblemBase` to write a custom C++
    class that provides the necessary methods, most importantly a
    `value` method that evaluates your function, and a `gradient`
    method that evaluates its gradient. See
    `dune/nonlinopt/userinterface.hh` for the definition of the
    base class.
* A `FiniteDifferenceProblemBase` class is available in the same
    file if you would like to avoid implementing your own gradients.
    This class uses central finite differences on your value method
    to compute gradients. Note that this is typically very
    inefficient, so define your gradients in a different way if that
    is possible.
* The two Boolean `subsequent` flags of the methods and the combined
    `value_and_grad` method can be used to reuse intermediate
    results and achieve a more efficient problem definition. See the
    documentation of these methods for details.
* The `dim` method is only used to detect vectors that are
    default-constructed and need to be initialized by the `zero` method,
    the solver itself is dimension-agnostic. Simply return a non-zero
    constant if you have a special case where the dimension of your
    problem is not really a well-defined concept.
* You can provide simple defaults for the `zero` and `hook` methods
    if you have no need for anything special, see the `ParsedProblem`
    class for suggestions.


Choosing an Optimization Method
-------------------------------

* The default configuration of the solver class provides an L-BFGS
    method with a window size of ten and directional scaling enabled.
    This algorithm performs reasonably well on a large range of
    problems, so try this one first. You can control the amount of
    program output using the `verbosity` constructor argument.
* Call `set_cg()` on the solver to configure it as a nonlinear
    conjugate gradients (NCG) method instead. This will produce the
    Hager-Zhang variant of nonlinear CG. See below on how to switch
    to a different version.
* Call `set_bfgs_cg()` to configure the solver as a variable-metric
    preconditioned nonlinear CG method, using the L-BFGS matrix as
    preconditioner (SCG). By default, this is a Hestenes-Stiefel CG
    method. Again, see below on how to change that if needed.
* Call `set_gmres()` on the solver to configure it as a nonlinear
    generalized minimal residual (N-GMRES) method. This method uses
    another solver, steepest descent by default, as a `preconditioner`,
    and extrapolates from its iterations based on a local least
    squares problem. There are also other variants of N-GMRES available.
* You can mix-and-match different configuration options by calling
    the `set_xyz<ABC>(...)` methods, where `xyz` is the part of the
    method you would like to replace (e.g., line search, termination
    criterion, type of nonlinear CG), and `ABC` is the class you
    would like to use. This class can be user-supplied, i.e., you
    can provide your own custom algorithms and criteria. The
    arguments of the methods are forwarded to the constructor call
    of the relevant class.


Optional: Provide a Preconditioner
----------------------------------

* You can provide a preconditioner using the `set_preconditioner(...)`
    method. This preconditioner is then employed in each iteration,
    and the preconditioned gradient is used instead of the gradient
    when determining the step direction.


Solving the Optimization Problem
--------------------------------

* You can call the `report()` method to get a print-out of the
    current high-level configuration and check that it meets your
    expectations.
* Simply call the `apply(problem,point)` method with an instance of
    your problem definition and the desired starting position. The
    solver will try to find a local minimum, and return it through
    its argument.
* The method returns the local minimizer it found through modification
    of the function argument `point`. There are `get_value` and
    `get_gradient` methods available if you need the final function
    value and gradient, with a Boolean flag for the latter if you are
    interested in the preconditioned gradient.
* If you need more control over the process, e.g., when using
    dune-nonlinopt within a larger simulation framework, you may
    call `hard_reset(problem)` to set up the solver without directly
    starting the optimization run. You can then call the `step()`
    method directly to perform a single step of the method, perform
    any other tasks, and call `step()` again when you are ready.


CUTEst Integration
==================

dune-nonlinopt can be used to run the [CUTEst benchmark][13] (or,
to be more precise, the unconstrained part thereof). See the file
[cutest/README.md](./cutest/README.md) for details.


References
==========

* Liu, D.C., and Nocedal, J.: On the limited memory BFGS method for
    large scale optimization, [(link)][1]
* De Sterck, H.: A Nonlinear GMRES Optimization Algorithm for
    Canonical Tensor Decomposition, [(link)][2]
* Hager, W.W., and Zhang, H.: A Survey of Nonlinear Conjugate
    Gradient Methods, [(link)][3]
* Barzilai, J., and Borwein, J.M.: Two-Point Step Size Gradient
    Methods, [(link)][4]
* Nocedal, J.: Updating quasi-Newton matrices with limited
    storage, [(link)][5]
* Hager, W.W., and Zhang, H.: A Nonmonotone Line Search Technique
    and Its Application to Unconstrained Optimization, [(link)][6]
* Hager, W.W., and Zhang, H.: Algorithm 851: CG_DESCENT, a conjugate
    gradient method with guaranteed descent, [(link)][7]
* Moré, J.J., Garbow, B.S., and Hillstrom, K.E.: Testing
    Unconstrained Optimization Software, [(link)][8]


[1]:  https://doi.org/10.1007/BF01589116
[2]:  https://arxiv.org/abs/1105.5331
[3]:  http://www.caam.rice.edu/~yzhang/caam554/pdf/cgsurvey.pdf
[4]:  https://doi.org/10.1093/imanum/8.1.141
[5]:  https://doi.org/10.1090/S0025-5718-1980-0572855-7
[6]:  https://doi.org/10.1137/S1052623403428208
[7]:  https://doi.org/10.1145/1132973.1132979
[8]:  https://dune-project.org
[9]:  https://eigen.tuxfamily.org/index.php?title=Main_Page
[10]: https://gitlab.dune-project.org/core/dune-common
[11]: https://github.com/beltoforion/muparser
[12]: https://doi.org/10.1145/355934.355936
[13]: https://www.cuter.rl.ac.uk/Problems/mastsif.shtml
