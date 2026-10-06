#!/bin/bash
# Installation choices for the comparison's private .local environment.
# Paths are relative to this repository. Existing code checkouts are reused;
# setup checks their revisions without fetching or changing their branches.

export PYTHON_VERSION="3.11"
export COCOA_PATH="../cocoa/Cocoa"
export ONECOV_PATH="../OneCovariance"
export ONECOV_GIT_COMMIT="311c2cfbe9d584d9d29abd3d8edaeea673470b7f"
export CAMB_GIT_COMMIT="6c978c716e163157c95191ee1c089e8b36c8adb9"
export HEALPY_VERSION="1.19.0"

# Match the pilot's numerical packages. --no-dependencies prevents pip
# from replacing the Conda base or silently selecting a newer NumPy.
PIPCP=(
  "numpy==1.26.3"
  "scipy==1.12.0"
  "astropy==6.1.7"
  "matplotlib==3.10.1"
  "hmf==3.5.2"
  "pyyaml==6.0.2"
  "pybind11==2.13.6"
  "setuptools==80.3.1"
  "setuptools-scm==8.3.1"
  "wheel==0.45.1"
  "cython==3.0.12"
  "pykg-config==1.3.0"
  "packaging==25.0"
  "sympy==1.14.0"
  "mpmath==1.3.0"
  "click==8.1.8"
  "toml==0.10.2"
  "deprecation==2.1.0"
  "cached_property==2.0.1"
  "rich==15.0.0"
  "markdown-it-py==4.2.0"
  "mdurl==0.1.2"
  "pygments==2.19.1"
)
