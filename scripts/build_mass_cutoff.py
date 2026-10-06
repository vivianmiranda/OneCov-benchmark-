"""Build an isolated LSST interface with a selected halo-table mass domain.

The lower limit changes in a private copy of structs.c. An optional private
cosmo3D.c copy extends the sigma-integral k cutoff for sub-solar halo tests.
Production sources, installed interface, fits and correction stay untouched.
Run after sourcing start_cocoa.sh; compilation uses the project Makefile.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

from common import revision, sha256
from split_mass_rule import private_source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--log10-min", type=int, default=2)
    parser.add_argument("--sigma-log10-kmax", type=int, choices=(5, 7, 9, 15),
                        default=5)
    parser.add_argument("--split-tail", action="store_true")
    args = parser.parse_args()
    if args.output.exists() or args.log10_min not in (-20, -3, 2, 4, 6):
        parser.error("use a new directory and log10-min -20, -3, 2, 4 or 6")
    root = Path(os.environ["ROOTDIR"])
    core = root/"external_modules/code/cosmolike_core"
    project = root/"projects/lsst_y1/interface"
    target = args.output.resolve()
    target.mkdir(parents=True)

    # Retain every production algorithm. This initializer controls the
    # supported sigma-table domain, not the chosen mass-integration edges.
    source = core/"cosmolike/structs.c"
    before = ".halo_m = {1.0e+4, 1.0e+17}"
    after = f".halo_m = {{1.0e{args.log10_min:+d}, 1.0e+17}}"
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
    sigma_source = core/"cosmolike/cosmo3D.c"
    sigma_record = dict(kmax_hmpc=1e5, source_sha256=sha256(sigma_source))
    if args.sigma_log10_kmax != 5:
        # The power reader still extrapolates the same supplied spectrum.
        # Changing its integration bound tests numerical tail coverage,
        # not the physical validity of that small-scale continuation.
        sigma_before = "const double kmax = 1.e5;"
        sigma_after = f"const double kmax = 1.e{args.sigma_log10_kmax};"
        sigma_content = sigma_source.read_text()
        sigma_original = "${ROOTDIR}/external_modules/code/cosmolike/cosmo3D.c"
        if sigma_content.count(sigma_before) != 1:
            raise ValueError("the production sigma cutoff has changed")
        if makefile.count(sigma_original) != 1:
            raise ValueError("the production sigma source list has changed")
        (target/"cosmo3D.c").write_text(
            sigma_content.replace(sigma_before, sigma_after))
        makefile = makefile.replace(sigma_original, "./cosmo3D.c")
        sigma_record.update(
            kmax_hmpc=10.0**args.sigma_log10_kmax,
            diagnostic_sha256=sha256(target/"cosmo3D.c"),
            replacement=dict(before=sigma_before, after=sigma_after),
        )
    makefile += f"\nCFLAGS += -I {core}/cosmolike\n"
    split_record = None
    if args.split_tail:
        halo_source = core/"cosmolike/covariances/halo_cov.c"
        (target/"halo_cov.c").write_text(private_source(halo_source.read_text()))
        original = "${ROOTDIR}/external_modules/code/cosmolike/covariances/halo_cov.c"
        if makefile.count(original) != 1:
            raise ValueError("the native halo source list has changed")
        makefile = makefile.replace(original, "./halo_cov.c")
        makefile += f"CFLAGS += -I {core}/cosmolike/covariances\n"
        split_record = dict(
            source_sha256=sha256(halo_source),
            diagnostic_sha256=sha256(target/"halo_cov.c"),
            rule_environment="COCOA_DIAGNOSTIC_TAIL_NQUAD",
            split_mass=1e4,
        )
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
        sigma_integral=sigma_record,
        split_tail=split_record,
        interface_sha256=sha256(target/"cosmolike_lsst_y1_interface.so"),
        production_interface_sha256=sha256(
            project/"cosmolike_lsst_y1_interface.so"),
    )
    (target/"build.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"Built isolated mass-domain diagnostic: {target}")


if __name__ == "__main__":
    main()
