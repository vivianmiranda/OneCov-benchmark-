#!/bin/bash
# Prepare the comparison environment, using Cocoa's Conda base.
# Source from this repository in a fresh Bash terminal after activating
# the Cocoa Conda environment. Downloads and pip setup happen here;
# compile_onecov.sh builds the native extensions without network access.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source setup_onecov.sh" >&2
  exit 1
fi

# Keep installation variables and directory changes out of the caller's
# shell. This does not activate .local; start_onecov.sh does that later.
(
  cd "$(dirname "${BASH_SOURCE[0]}")" || return 1
  source ./set_installation_options.sh || return 1

  if [[ -z "${CONDA_PREFIX:-}" || -n "${VIRTUAL_ENV:-}" ||
        -n "${ROOTDIR:-}" ]]; then
    echo "Use a fresh Bash terminal with only 'conda activate cocoa'." >&2
    return 1
  fi

  PYTHON="${CONDA_PREFIX:?}/bin/python"
  if [[ "$("${PYTHON}" -c \
    'import sys; print("%d.%d" % sys.version_info[:2])')" != \
    "${PYTHON_VERSION}" ]]; then
    echo "The Cocoa Conda base must use Python ${PYTHON_VERSION}." >&2
    return 1
  fi

  # The comparison uses the installed Cocoa CAMB, including its patches.
  # Reusing that source keeps the Boltzmann code common to both sides.
  CAMB_PATH="${COCOA_PATH}/external_modules/code/CAMB"
  if [[ "$(git -C "${ONECOV_PATH}" rev-parse HEAD)" != \
        "${ONECOV_GIT_COMMIT}" ||
        "$(git -C "${CAMB_PATH}" rev-parse HEAD)" != \
        "${CAMB_GIT_COMMIT}" ]]; then
    echo "Checkout revisions differ from set_installation_options.sh." >&2
    return 1
  fi
  if [[ ! -f "${CAMB_PATH}/camb/camblib.so" ]]; then
    echo "Compile CAMB through Cocoa before setting up the comparison." >&2
    return 1
  fi

  "${PYTHON}" -m venv --system-site-packages .local || return 1
  PYTHON="$(pwd -P)/.local/bin/python"

  # Pin installations inside .local. Packages already supplied at the
  # requested versions by the Conda base are reused, never upgraded there.
  "${PYTHON}" -m pip install --no-dependencies \
    --prefer-binary "${PIPCP[@]}" || return 1

  "${PYTHON}" - "${CAMB_PATH}" <<'PY' || return 1
import sys
import sysconfig
from pathlib import Path

camb = Path(sys.argv[1]).resolve()
site = Path(sysconfig.get_path("purelib"))
(site / "benchmark_camb.pth").write_text(
    f"import sys; sys.path.insert(0, {str(camb)!r})\n"
)
PY

  # Keep healpy's build tree at a permanent location: on macOS its
  # extension can reference libraries built inside that tree.
  mkdir -p external_modules/code external_modules/data || return 1
  ARCHIVE="external_modules/data/healpy-${HEALPY_VERSION}.tar.gz"
  if [[ ! -f "${ARCHIVE}" ]]; then
    "${PYTHON}" -m pip download --no-dependencies --no-build-isolation \
      --no-binary=:all: --dest external_modules/data \
      "healpy==${HEALPY_VERSION}" || return 1
  fi
  if [[ ! -d "external_modules/code/healpy-${HEALPY_VERSION}" ]]; then
    tar -xzf "${ARCHIVE}" -C external_modules/code || return 1
  fi

  echo "Setup complete. Next: source compile_onecov.sh"
) || return 1

return 0
