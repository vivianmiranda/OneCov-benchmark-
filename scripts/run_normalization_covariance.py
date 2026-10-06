"""Test full LSST Y1 covariance sensitivity to halo normalization.

Run after start_cocoa.sh. Choose native, mass_only or mass_and_bias. The
last two use an explicit supplied-moment adapter; neither installed source
nor the fitted functions exposed by the ordinary interface are changed.
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
from halo_normalization import NormalizedMoments, check_ingredients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("native", "mass_only", "mass_and_bias",
                                         "ingredients"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error("choose a new output directory")
    signal.alarm(900)

    root = Path(os.environ["ROOTDIR"])
    core = root/"external_modules/code/cosmolike_core"
    project = root/"projects/lsst_y1"
    sys.path[:0] = [str(project/"interface"), str(core),
                   str(project/"covariance")]
    import cosmolike_lsst_y1_interface as ci
    import lsst_y1_covariance as survey
    from cosmolike_notebook_utils.covariance.command_line import (
        load_run_configuration,
    )
    from cosmolike_notebook_utils.covariance.forecast import save_forecast

    os.chdir(root)
    filename = project/"EXAMPLE_EVALUATE_COVARIANCE.yaml"
    settings, run = load_run_configuration(filename, survey)
    settings["execution"] = dict(
        backend="production with explicit diagnostic moment adapter",
        threads=run["threads"],
        halo_normalization=args.mode,
    )
    if run["camb_path"] is not None:
        sys.path.insert(0, str(Path(run["camb_path"]).resolve()))
    ci.set_omp_threads(n=run["threads"])
    output.mkdir(parents=True)
    start = time.perf_counter()
    tables = survey.initialize(interface=ci, settings=settings)
    initialization = time.perf_counter()-start
    ci.set_omp_threads(n=run["threads"])

    if args.mode == "ingredients":
        record = check_ingredients(
            interface=ci, edges=settings["lnm_edges"],
            cosmology=settings["cosmology"],
        )
        record["core"] = revision(core)
        record["interface_sha256"] = sha256(ci.__file__)
        record["adapter_sha256"] = sha256(
            Path(__file__).with_name("halo_normalization.py"))
        (output/"report.json").write_text(json.dumps(record, indent=2)+"\n")
        print("Normalization and direct mass-integral checks passed.")
        return

    backend = ci.covariance
    if args.mode != "native":
        backend = NormalizedMoments(interface=ci, mode=args.mode)

    def progress(stage, seconds):
        print(f"{stage}: {seconds:.3f} s", flush=True)

    result = survey.compute(
        interface=ci, settings=settings, space=run["space"],
        progress=progress, backend=backend,
    )
    save_forecast(result, output/"covariance.npz")
    np.savez_compressed(output/"power_tables.npz", **tables)
    record = dict(
        mode=args.mode,
        core=revision(core),
        project=revision(project),
        interface_sha256=sha256(ci.__file__),
        yaml_sha256=sha256(filename),
        mass_min=float(np.exp(settings["lnm_edges"][0])),
        shape=list(result["total"].shape),
        stages_s=result["stages_s"],
        initialization_s=initialization,
        threads=run["threads"],
        covariance_sha256=sha256(output/"covariance.npz"),
        power_table_sha256=sha256(output/"power_tables.npz"),
        script_sha256=sha256(__file__),
        adapter_sha256=sha256(Path(__file__).with_name("halo_normalization.py")),
    )
    (output/"report.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"Saved complete covariance: {output}", flush=True)


if __name__ == "__main__":
    main()
