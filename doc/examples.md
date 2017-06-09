## Examples

To run the examples, cuIBM should be compiled with the executable `cuibm` available in the `bin` folder of the software directory `$CUIBM_DIR` (see [installation instructions](installation.md)).
You can add the binary executable to your `PATH` environment variable:

    > export PATH=$CUIBM_DIR/bin:$PATH

Examples:

* [lid-driven cavity flow](lidDrivenCavity.md) at Reynolds numbers 100, 1000, 3200, and 5000;
* [impulsively started cylinder](cylinder.md) at Reynolds numbers 40, 550, and 3000;
* [flapping wing](flapping.md) at Reynolds number 75.

Python post-processing scripts are available each case folder; they make use of the package [`snake`](https://github.com/mesnardo/snake) that is bundled in the `external` folder of cuIBM (version `0.3`).
To install `snake`:

    > cd $CUIBM_DIR/external/snake-0.3
    > python setup.py install

To post-process the numerical solution from cuIBM, `snake` requires Python (2.7 or 3.5), Matplotlib, Scipy, and Pandas.
