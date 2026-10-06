#!/bin/bash
# Deactivate only this comparison environment; retain the Conda base.
# The virtual environment's deactivate restores PATH and the prompt.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source stop_onecov.sh" >&2
  exit 1
fi

if [[ "${VIRTUAL_ENV:-}" != \
      "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/.local" ]]; then
  echo "This comparison's private environment is not active." >&2
  return 1
fi

deactivate || return 1
return 0
