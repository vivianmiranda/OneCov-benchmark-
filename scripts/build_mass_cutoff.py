"""Build an isolated LSST interface with a wider halo-table mass domain.

Only the lower limit in a private copy of structs.c changes. The production
sources, installed interface, fitted functions and correction stay untouched.
Run after sourcing start_cocoa.sh; compilation uses the project Makefile.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

from common import revision, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--log10-min", type=int, default=2)
    args = parser.parse_args()
    if args.output.exists() or args.log10_min not in (2, 4, 6):
        parser.error("use a new output directory and log10-min of 2, 4 or 6")
    root = Path(os.environ["ROOTDIR"])
    core = root/"external_modules/code/cosmolike_core"
    project = root/"projects/lsst_y1/interface"
    target = args.output.resolve()
    target.mkdir(parents=True)

    # Retain every production algorithm. This initializer controls the
    # supported sigma-table domain, not the chosen mass-integration edges.
    source = core/"cosmolike/structs.c"
    before = ".halo_m = {1.0e+6, 1.0e+17}"
    after = f".halo_m = {{1.0e+{args.log10_min}, 1.0e+17}}"
    content = source.read_text()
    if content.count(before) != 1:
        raise ValueError("the production mass initializer has changed")
    (target/"structs.c").write_text(content.replace(before, after))
    shutil.copy2(project/"interface.cpp", target/"interface.cpp")

    makefile = (project/"MakefileCosmolike").read_text()
    original = "${ROOTDIR}/external_modules/code/cosmolike/structs.c"
    if makefile.count(original) != 1:
        raise ValueError("the production build source list has changed")
    makefile = makefile.replace(original, "./structs.c")
    makefile += f"\nCFLAGS += -I {core}/cosmolike\n"
    (target/"MakefileCosmolike").write_text(makefile)
    environment = dict(os.environ)
    environment.pop("IGNORE_COSMOLIKE_LSST_Y1_COVARIANCE", None)
    environment.pop("COSMOLIKE_DEBUG_MODE", None)
    environment.pop("COSMOLIKE_AGGRESSIVE_MODE", None)
    with (target/"build.log").open("w") as log:
        subprocess.run(
            ["make", "-j1", "-f", "MakefileCosmolike", "all"],
            cwd=target, env=environment, stdout=log,
            stderr=subprocess.STDOUT, check=True, timeout=600,
        )
    record = dict(
        core=revision(core), project=revision(project),
        table_mass_min=10.0**args.log10_min, table_mass_max=1e17,
        production_structs_sha256=sha256(source),
        diagnostic_structs_sha256=sha256(target/"structs.c"),
        replacement=dict(before=before, after=after),
        interface_sha256=sha256(target/"cosmolike_lsst_y1_interface.so"),
        production_interface_sha256=sha256(
            project/"cosmolike_lsst_y1_interface.so"),
    )
    (target/"build.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"Built isolated mass-domain diagnostic: {target}")


if __name__ == "__main__":
    main()
