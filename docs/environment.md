# Environments for the comparison

The input exporter runs in Cocoa's existing environment, with its LSST Y1
covariance bindings enabled. The OneCovariance runner can use a separate
environment. Files in `work/` connect the two; neither environment needs
to upgrade the other.

## Dependencies

Use Python 3.11 for the currently tested setup. OneCovariance requires
NumPy, SciPy, astropy, hmf, CAMB, healpy and matplotlib, plus its **bundled**
Levin extension. GSL and a working C++/OpenMP toolchain are needed to build
that extension. An unrelated package with the same `levin` import name
is not a substitute.

Follow the installation instructions in the local OneCovariance checkout
to prepare an independent environment. To build only its bundled extension
once dependencies are present, run from this benchmark root:

```bash
python -m pip install --no-deps --no-build-isolation ../OneCovariance
```

This installs the extension; the runner still needs `--onecov` to locate
the source checkout. Test imports before starting a covariance:

```bash
python -c "import camb, healpy, hmf, levin; print(camb.__version__); print(levin.Levin)"
```

OneCovariance imports several estimator modules even for a Gaussian
Fourier run. A missing healpy dependency can therefore appear at the end
of a traceback as a misleading missing `cov_theta_space` import. Inspect
the **first** exception in the traceback.

## Setup used for the initial tests

The isolated `.venv` inherits the existing Python 3.11 Cocoa environment's
NumPy/SciPy/astropy packages without modifying them. New packages were
installed only in `.venv`. The numerical versions are:

| Package | Version |
| --- | --- |
| NumPy | 1.26.3 |
| SciPy | 1.12.0 |
| astropy | 6.1.7 |
| hmf | 3.5.2 |
| healpy | 1.19.0 |
| CAMB | 1.6.7, from Cocoa's local source checkout |
| Levin | 0.0.1, built from OneCovariance's bundled source |

The CAMB version is an import from the same source as Cocoa, not a PyPI
release pin. For an isolated environment that should use an existing CAMB
checkout, put its parent directory on `PYTHONPATH` before the run, or add
that directory to a `.pth` file in the isolated environment. Do not copy
CAMB arrays from a different cosmology to compensate for a version mismatch.
`run.json` records the actual CAMB version and import path after preflight.

On this macOS installation, the first healpy source wheel referenced
libraries in a removed temporary build directory. Rebuilding from an
unmodified source archive in a persistent directory resolved the imports.
Keep that build directory while the installed extension depends on it.
This is an environment issue; no OneCovariance numerical source was changed.

## Running sequentially

After activating the chosen environment, set the worker count explicitly:

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=disabled
```

The runner derives OneCovariance's `num_cores` from `OMP_NUM_THREADS` and
fixes BLAS to one thread before imports. There is no second thread count
in the survey configuration. Use one numerical job at a time; eight
workers are the laptop limit, not a claim that every OneCov stage uses
all eight.

Logs and a wall-clock deadline are safeguards for affordable accuracy
tests. Initial feasibility times are not a performance comparison. Do not
optimize OneCovariance in response to an expensive stage; reduce the
number of compared bins or ingredients while preserving the accuracy
needed for that particular test.
