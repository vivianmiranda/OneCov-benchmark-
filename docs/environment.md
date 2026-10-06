# Environments for the comparison

The input exporter runs in Cocoa's environment, with its LSST Y1 covariance
bindings enabled. OneCovariance runs in this repository's `.local`
environment. Both use the Cocoa Conda base; files in `work/` connect their
calculations. Open separate Bash terminals for the two runtime environments.

For the Cocoa export, use the activation, covariance compilation and
platform-specific thread settings in [Step 1 of the README](../README.md#reproduction).
These follow the main Cocoa and LSST Y1 READMEs; no additional Cocoa
environment variables or dependency installation are required here.

## Installation choices

Use Python 3.11 for the currently tested setup. OneCovariance requires
NumPy, SciPy, astropy, hmf, CAMB, healpy and matplotlib, plus its **bundled**
Levin extension. GSL and a working C++/OpenMP toolchain are needed to build
that extension. An unrelated package with the same `levin` import name
is not a substitute.

Use the [installation steps](../README.md#installation). Choices belong to
`set_installation_options.sh`, following Cocoa's setup convention.

| Choice | Meaning |
| --- | --- |
| `COCOA_PATH` | Installed Cocoa runtime directory containing CAMB. |
| `ONECOV_PATH` | Existing OneCovariance source checkout. |
| `CAMB_GIT_COMMIT`, `ONECOV_GIT_COMMIT` | Revisions checked before setup. |
| `HEALPY_VERSION` | Source version downloaded during setup. |
| `PIPCP` | Pinned Python packages installed only in `.local`. |

`setup_onecov.sh` may download packages. It reuses the Conda base through
Python's `--system-site-packages` option and installs missing or different
pinned versions inside `.local`. It neither changes the Conda base nor
uses Cocoa's private `.local` directory.

`compile_onecov.sh` builds the two extensions with `--no-index`,
`--no-dependencies` and `--no-build-isolation`. It uses the compilers and
GSL from the activated Conda environment. Successful compilation ends with
an import check of CAMB, healpy, hmf, Levin and the OneCovariance reader.

The healpy source and build tree remain in `external_modules/code` because
its macOS extension can depend on libraries in that tree. Keep this
directory alongside `.local`. Both directories are ignored by Git.

> [!TIP]
> If an import fails with `cov_theta_space` missing, inspect the first
> exception: OneCovariance loads several estimator modules even for a
> Gaussian Fourier case, and a missing healpy library can cause that error.

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

The setup script registers Cocoa's compiled CAMB source in `.local`, so
both sides use the same Boltzmann code and its Cocoa installation patches.
It does not install a different CAMB release from PyPI. `run.json` records
the actual CAMB version and import path after preflight.

The `.venv` mentioned in the initial result record is retained as pilot
provenance. The documented setup scripts create `.local` independently;
they do not rename or overwrite that environment or the saved results.

## Running sequentially

After `conda activate cocoa` and `source start_onecov.sh`, set the worker
count explicitly:

```bash
export OMP_NUM_THREADS=8
```

Keep the platform-specific OpenMP settings from Step 1 of the README:
`OMP_PROC_BIND=close` on Linux or `disabled` on macOS (arm), with
`OMP_PLACES=cores` and `OMP_DYNAMIC=FALSE`.

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
