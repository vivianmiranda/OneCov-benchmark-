#!/bin/bash
# Build healpy and OneCovariance's bundled Levin extension into .local.
# Source after setup_onecov.sh, with the Cocoa Conda base active.
# This step uses local sources only and must not access the network.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source compile_onecov.sh" >&2
  exit 1
fi

(
  cd "$(dirname "${BASH_SOURCE[0]}")" || return 1
  source ./set_installation_options.sh || return 1

  if [[ -z "${CONDA_PREFIX:-}" || ! -x .local/bin/python ||
        -n "${VIRTUAL_ENV:-}" || -n "${ROOTDIR:-}" ]]; then
    echo "Activate only the Cocoa Conda base and run setup_onecov.sh first." >&2
    return 1
  fi

  # Conda supplies the compilers, GSL and OpenMP libraries. The runtime
  # search path allows the installed extensions to find those libraries.
  export PATH="$(pwd -P)/.local/bin:${PATH}"
  export CPPFLAGS="-I${CONDA_PREFIX}/include ${CPPFLAGS:-}"
  export LDFLAGS="-L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib ${LDFLAGS:-}"
  # Let healpy find CFITSIO's matching header and library in Conda.
  # Otherwise it can build a second copy but link to Conda's older copy.
  export PKG_CONFIG_PATH="${CONDA_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  PYTHON="$(pwd -P)/.local/bin/python"

  "${PYTHON}" -m pip install \
    "./external_modules/code/healpy-${HEALPY_VERSION}" \
    --no-dependencies --no-index --no-build-isolation || return 1

  "${PYTHON}" -m pip install "${ONECOV_PATH}" \
    --no-dependencies --no-index --no-build-isolation || return 1

  # Imports load the actual native libraries, so a missing GSL/OpenMP
  # library fails here before a user starts a covariance calculation.
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "${PYTHON}" - "${ONECOV_PATH}" <<'PY' || return 1
import sys
from pathlib import Path

import camb
import healpy
import hmf
import levin

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from onecov.cov_input import Input

assert hasattr(levin, "Levin")
print("CAMB:", camb.__version__, camb.__file__)
print("healpy:", healpy.__version__)
print("hmf:", hmf.__version__)
print("OneCovariance imports passed.")
PY

  echo "Compilation complete. Next: source start_onecov.sh"
) || return 1

return 0
