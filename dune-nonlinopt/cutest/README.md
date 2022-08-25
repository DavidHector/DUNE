CUTEst interface for dune-nonlinopt
===================================

These subfolders contain the interface code that is needed
to use dune-nonlinopt together with the CUTEst Test Problem
Set. If you would like to run CUTEst with dune-nonlinopt,
proceed as follows:

* Install the lastest [CUTEst test suite][1], following the
    official instructions.
* You should now have a folder named `cutest` somewhere, with
    subfolders `packages` and `src`. These folders are mirrored
    in the current directory. Place the file under `packages/defaults`
    in the CUTEst folder with the same name, this makes
    dune-nonlinopt known as a package in CUTEst.
* Copy the folder `src/dune-nonlinopt` into the `src` folder
    of CUTEst. This folder contains the interface code, which
    defines a problem class wrapper and a simple main function
    that calls CUTEst, and a corresponding makefile.
* Modify the variable `NONLINOPT_PATH` in the file `makemaster`,
    and point it to the directory containing the header files.
    CUTEst should now be able to build the interface code, link
    it against the FORTRAN code providing problem definitions,
    and run the resulting program.
* You can run different methods by modifying the interface
    code and using the appropriate function calls as explained
    in the main README.

[1]: https://www.cuter.rl.ac.uk/Problems/mastsif.shtml
