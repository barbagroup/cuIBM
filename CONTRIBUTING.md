# How to contribute to cuIBM

Welcome to the developer's guide of cuIBM!

## Adding new features and fixing bugs

All new features and bug fixes must go through a pull-request review procedure.
If you want to add something to cuIBM, please fork the Barbagroup's [cuIBM](https://github.com/barbagroup/cuIBM) repository, make your changes on your fork, and then open a pull-request with your changes against the main cuIBM repository.

For new features and minor bugs (with small impact), the base branch of the pull-request should be the `develop` branch of the main repository.
(The `develop` branch will be merged into the `master` one once we are ready for a new release of cuIBM.)

For major bugs, the base branch should be the `master` one of the main repository; it will be considered as a hotfix (bugfix) and a new version of cuIBM will be released as soon as possible by the maintainers with the micro number incremented.

New features should come with some kind of test or example to verify and/or validate the implementation.
For example, the Markdown file [`doc/flapping.md`](https://github.com/barbagroup/cuIBM/blob/master/doc/flapping.md) provides a detailed description of a flapping-wing simulation at Reynolds number 75; it serves as validation for the case of a single moving body (using the immersed-boundary projection method) and compares cuIBM results with ones from other studies (experimental and computational).


## Reporting bugs and requesting new features

To report bugs, request new features, or simply ask questions, please open a GitHub issue on the Barbagroup's cuIBM repository.


## Writing documentation

New classes, methods, functions, and namespaces must be documented with Doxygen-style doctrings.
(See [here](https://github.com/barbagroup/cuIBM/blob/master/src/solvers/NavierStokesSolver.h) for reference.)

To locally check the Doxygen documentation, use the command-line:

    > make doc

and open the file `doc/html/index.html` in your favorite browser.

You should also add code documentation whenever necessary; it will greatly help other developers to review your new features and bug fixes.

For new features, user's documentation must also be written.
For this purpose, we use Markdown files that are located in the `doc` folder of the root directory of cuIBM.

The Wiki pages and Doxygen API documentation are up-to-date with the latest release of cuIBM, which should also be the latest commit on the `master` branch.

Once a new new release is drafted, we will merge the `doc` folder into the Wiki of cuIBM and the API documentation will be merged into the branch `gh-pages`.
