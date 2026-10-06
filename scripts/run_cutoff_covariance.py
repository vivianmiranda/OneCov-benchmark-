"""Compute a complete LSST Y1 covariance for one halo-mass cutoff.

The official YAML reader, survey setup and production covariance backend
are reused. Only the mass integration edges change. An optional isolated
interface supplies the wider sigma-table domain; no installed files change.
Run after sourcing start_cocoa.sh, with OMP_NUM_THREADS set in the shell.
"""

import argparse
import json
import os
from pathlib import Path
import signal
import sys
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["COBAYA_NOMPI"] = "1"

import numpy as np

from common import revision, sha256
from diagnose_mass_cutoff import mass_edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interface", type=Path)
    parser.add_argument("--log10-min", type=int, choices=(2, 4, 6), default=6)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error("choose a new output directory")
    if args.interface is None and args.log10_min != 6:
        parser.error("lower cutoffs require the isolated wider-domain build")
    selected = None if args.interface is None else args.interface.resolve()
    output.mkdir(parents=True)
    signal.alarm(900)

    root = Path(os.environ["ROOTDIR"])
    core = root/"external_modules/code/cosmolike_core"
    project = root/"projects/lsst_y1"
    location = project/"interface" if selected is None else selected
    sys.path[:0] = [str(location), str(core), str(project/"covariance")]
    import cosmolike_lsst_y1_interface as ci
    import lsst_y1_covariance as survey
    from cosmolike_notebook_utils.covariance.command_line import (
        load_run_configuration,
    )
    from cosmolike_notebook_utils.covariance.forecast import save_forecast

    build = None
    if selected is not None:
        build = json.loads((selected/"build.json").read_text())
        if sha256(ci.__file__) != build["interface_sha256"]:
            raise ValueError("imported interface differs from diagnostic build")

    # The example's relative theory paths follow the ordinary CLI contract.
    # Preserve its Gaussian non-Limber choice and every accuracy setting.
    os.chdir(root)
    filename = project/"EXAMPLE_EVALUATE_COVARIANCE.yaml"
    settings, run = load_run_configuration(filename, survey)
    settings["lnm_edges"] = mass_edges(args.log10_min)
    settings["execution"] = dict(
        backend="production", threads=run["threads"],
        diagnostic="mass cutoff with unchanged halo prescriptions",
        table_domain="native" if build is None else "extended to 1e2",
    )
    if run["camb_path"] is not None:
        sys.path.insert(0, str(Path(run["camb_path"]).resolve()))
    ci.set_omp_threads(n=run["threads"])
    start = time.perf_counter()
    tables = survey.initialize(interface=ci, settings=settings)
    initialization = time.perf_counter()-start
    ci.set_omp_threads(n=run["threads"])
    print(f"Initialization including CAMB: {initialization:.3f} s", flush=True)

    def progress(stage, seconds):
        print(f"{stage}: {seconds:.3f} s", flush=True)

    result = survey.compute(interface=ci, settings=settings, space=run["space"],
                            progress=progress, backend=ci.covariance)
    save_forecast(result, output/"covariance.npz")
    np.savez_compressed(output/"power_tables.npz", **tables)
    record = dict(
        core=revision(core), project=revision(project),
        yaml_sha256=sha256(filename), interface_sha256=sha256(ci.__file__),
        diagnostic_build=build, mass_min=10.0**args.log10_min,
        shape=list(result["total"].shape), stages_s=result["stages_s"],
        initialization_s=initialization, threads=run["threads"],
        covariance_sha256=sha256(output/"covariance.npz"),
        power_table_sha256=sha256(output/"power_tables.npz"),
        script_sha256=sha256(__file__),
    )
    (output/"report.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"Saved complete covariance: {output}", flush=True)


if __name__ == "__main__":
    main()
