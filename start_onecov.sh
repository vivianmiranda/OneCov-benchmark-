#!/bin/bash
# Activate the comparison's private Python environment in this shell.
# Like start_cocoa.sh, this is sourced after activating the Conda base.
# Use a separate terminal for Cocoa, so its runtime paths do not mix here.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source start_onecov.sh" >&2
  exit 1
fi

if [[ -z "${CONDA_PREFIX:-}" || -n "${VIRTUAL_ENV:-}" ||
      -n "${ROOTDIR:-}" ]]; then
  echo "Use a fresh Bash terminal with only 'conda activate cocoa'." >&2
  return 1
fi

source "$(dirname "${BASH_SOURCE[0]}")/.local/bin/activate" || return 1
return 0
