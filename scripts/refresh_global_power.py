"""Refresh small cross-code tests after global covariance power refinement.

Print the sequential commands by default. --execute runs them in new folders
and stops on the first failure. This is an accuracy refresh, not a timing
benchmark: initialization and numerical settings remain visible in child
records, but their elapsed times must not become new performance claims.

Start from the documented active Cocoa environment. The OneCov children use
its separate interpreter with Cocoa's Python import paths removed. No child
changes either library. Existing saved results are never overwritten.
"""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import time

from common import revision, sha256


def stages(args):
    """Build the existing component and complete-matrix comparison commands."""
    root = args.output
    inputs = root/"inputs"
    figures = root/"results/figures"
    seed = root/"native_seed/onecov.ini"
    commands = []

    def add(name, backend, script, *arguments, seconds=660):
        interpreter = getattr(args, f"{backend}_python")
        command = [str(interpreter), str(args.repo/"scripts"/script)]
        command.extend(str(value) for value in arguments)
        commands.append(dict(name=name, backend=backend, command=command,
                             deadline_seconds=seconds))

    add("inputs", "cocoa", "prepare_lsst_y1.py", "--cocoa", args.cocoa,
        "--output", inputs, "--integer-ell-range", 30, 1620,
        "--expected-power-nodes", 11993)

    # These Gaussian runs isolate assembly using the same angular spectra.
    # They cover shear, cross-probe blocks, and both low and high multipoles.
    for case, lower, upper in (("shear", 30, 150), ("3x2", 30, 150),
                               ("3x2", 1500, 1620)):
        name = f"{case}_{lower}_{upper}"
        add(name+"_onecov", "onecov", "run_onecov.py", "--onecov",
            args.onecov, "--inputs", inputs, "--output", root/name,
            "--case", case, "--terms", "gaussian", "--spectra",
            "shared-cells", "--band-limits", lower, upper, "--timeout", 600)
        add(name+"_compare", "cocoa", "compare_gaussian.py", root/name,
            "--cocoa", args.cocoa, "--output", root/"gaussian"/name)
    add("plot_gaussian", "cocoa", "plot_gaussian.py", "--comparisons",
        root/"gaussian", "--output", figures)

    # A fresh INI supplies survey and nonlinear-P inputs to native OneCov.
    # Its halo linear power still follows its own hmf/CAMB preparation.
    add("native_seed", "onecov", "run_onecov.py", "--onecov", args.onecov,
        "--inputs", inputs, "--output", seed.parent, "--case", "shear",
        "--spectra", "shared-power", "--prepare-only")
    add("ssc_response_export", "onecov", "ssc_response.py", "export", seed,
        "--onecov", args.onecov, "--output", root/"ssc_response_onecov")
    add("ssc_response_compare", "cocoa", "ssc_response.py", "compare",
        root/"ssc_response_onecov", "--output", root/"ssc_response_compare")
    for radial in (601, 1201):
        native = root/f"ssc_projection_{radial}"
        compared = root/f"ssc_comparison_{radial}"
        add(f"ssc_export_{radial}", "onecov", "compare_ssc.py", "export",
            seed, "--onecov", args.onecov, "--radial-nodes", radial,
            "--delta-z", .05, "--mass-nodes", 800, "--output", native,
            "--timeout", 600)
        add(f"ssc_compare_{radial}", "cocoa", "compare_ssc.py", "compare",
            native, "--cocoa", args.cocoa, "--output", compared)
    add("plot_ssc", "cocoa", "plot_ssc.py", root/"ssc_comparison_601",
        "--refinements", root/"ssc_comparison_1201", "--output", figures)

    # Keep AB=1 while varying only the quadrature. AB=2 also refines power,
    # so it would need its own shared input export, not this 11993-node one.
    for ingredient, script in (("halo", "compare_halo.py"),
                               ("trispectrum", "compare_trispectrum.py")):
        native = root/f"{ingredient}_onecov"
        add(ingredient+"_export", "onecov", script, "export", seed,
            "--onecov", args.onecov, "--mass-nodes", 800, "--output",
            native, "--timeout", 600)
        for level in (0, 2):
            add(f"{ingredient}_compare_i{level}", "cocoa", script,
                "compare", native, "--cocoa", args.cocoa,
                "--integration-accuracy", level, "--output",
                root/f"{ingredient}_cocoa_i{level}")
    add("plot_halo", "cocoa", "plot_halo.py", root/"halo_onecov",
        root/"halo_cocoa_i2", "--output", figures)
    add("connected_export", "onecov", "compare_connected.py", "export",
        seed, "--onecov", args.onecov, "--k-nodes", 129,
        "--mass-nodes", 400, "--radial-nodes", 601, "--delta-z", .125,
        "--output", root/"connected_onecov", "--timeout", 600)
    add("connected_compare", "cocoa", "compare_connected.py", "compare",
        root/"connected_onecov", "--output", root/"connected_cocoa")
    add("plot_connected", "cocoa", "plot_connected.py",
        root/"trispectrum_cocoa_i2", root/"connected_cocoa",
        "--output", figures)

    # Retain every G/SSC/cNG/total entry. Real space retains both xi signs
    # and their cross-covariance, using OneCov's shipped transform defaults.
    for space in ("fourier", "real"):
        for backend in ("onecov", "cocoa"):
            deadline = 1800 if space == "real" else 900
            add(f"{space}_{backend}", backend, "complete_shear.py", backend,
                "--space", space, "--inputs", inputs, "--cocoa", args.cocoa,
                "--onecov", args.onecov, "--output", root/f"{space}_{backend}",
                "--timeout", deadline, "--accuracy-only", seconds=deadline+60)
        add(f"plot_{space}", "cocoa", "plot_complete_shear.py",
            root/f"{space}_cocoa", root/f"{space}_onecov", "--output",
            root/f"results/{space}_shear.json", "--figure",
            figures/f"{space}_shear_difference.png")
    return commands


def check_power(root):
    """Require all completed native CoCoA tests to share installed power."""
    baseline = json.loads((root/"inputs/manifest.json").read_text())
    expected = baseline["installed_power"]
    if expected["nk"] != 11993 or expected["power_refinement"] != 8:
        raise ValueError("input export is not the adopted AB=1 power grid")
    reports = list(root.glob("*_cocoa_i*/comparison.json"))
    reports.extend(root.glob("*_cocoa/report.json"))
    for report in reports:
        current = json.loads(report.read_text())["installed_power"]
        if current != expected:
            raise ValueError(f"installed power differs from export: {report}")


def main():
    """Run a fresh, bounded sequence only when --execute is requested."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cocoa", type=Path, default=Path("../cocoa/Cocoa"))
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--cocoa-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--onecov-python", type=Path,
                        default=Path(".local/bin/python"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget-seconds", type=int, default=7200)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    args.repo = Path(__file__).resolve().parents[1]
    for key in ("cocoa", "onecov", "output"):
        setattr(args, key, getattr(args, key).resolve())
    # Keep venv interpreter symlinks: resolving them selects the base Python.
    for key in ("cocoa_python", "onecov_python"):
        setattr(args, key, getattr(args, key).absolute())
    commands = stages(args)
    if not args.execute:
        for stage in commands:
            print(stage["name"]+": "+shlex.join(stage["command"]))
        return
    if args.output.exists():
        parser.error("choose a new output folder; old results are preserved")
    if args.budget_seconds < 1:
        parser.error("the wall-time budget must be positive")
    workers = int(os.environ.get("OMP_NUM_THREADS", "0"))
    if not 1 <= workers <= 8 or os.environ.get("OPENBLAS_NUM_THREADS") != "1":
        parser.error("use 1–8 OMP_NUM_THREADS and OPENBLAS_NUM_THREADS=1")

    cocoa_env = os.environ.copy()
    onecov_env = cocoa_env.copy()
    for key in ("PYTHONPATH", "PYTHONHOME", "ROOTDIR", "VIRTUAL_ENV"):
        onecov_env.pop(key, None)
    core = args.cocoa/"external_modules/code/cosmolike_core"
    manifest = dict(
        status="running", purpose="accuracy only; not performance evidence",
        started_utc=datetime.now(timezone.utc).isoformat(),
        power_contract="initialize returns global natural-cubic refined tables",
        expected_power_nodes=11993, omp_threads=workers,
        budget_seconds=args.budget_seconds,
        revisions={"benchmark": revision(args.repo), "core": revision(core),
                   "onecov": revision(args.onecov)},
        scripts={p.name: sha256(p) for p in (args.repo/"scripts").glob("*.py")},
        stages=commands,
    )
    args.output.mkdir(parents=True)
    (args.output/"logs").mkdir()
    (args.output/"results/figures").mkdir(parents=True)
    ledger = args.output/"refresh.json"

    def save():
        ledger.write_text(json.dumps(manifest, indent=2)+"\n")

    started = time.monotonic()
    save()
    try:
        for stage in commands:
            remaining = args.budget_seconds-(time.monotonic()-started)
            if remaining < stage["deadline_seconds"]:
                raise TimeoutError("insufficient remaining budget for "
                                   +stage["name"])
            stage["status"] = "running"
            stage["log"] = f"logs/{stage['name']}.log"
            save()
            print(stage["name"], flush=True)
            environment = cocoa_env if stage["backend"] == "cocoa" else onecov_env
            with (args.output/stage["log"]).open("w") as log:
                process = subprocess.Popen(
                    stage["command"], cwd=args.repo, env=environment,
                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                )
                try:
                    code = process.wait(timeout=stage["deadline_seconds"])
                except BaseException:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait()
                    raise
            stage["returncode"] = code
            if code != 0:
                raise RuntimeError(f"{stage['name']} failed; see {stage['log']}")
            check_power(args.output)
            stage["status"] = "completed"
            save()
        manifest["status"] = "completed"
        manifest["reports"] = {
            str(path.relative_to(args.output)): sha256(path)
            for path in args.output.rglob("*.json") if path != ledger
        }
    except BaseException as error:
        manifest["status"] = "stopped"
        manifest["error"] = str(error)
        save()
        raise
    save()


if __name__ == "__main__":
    main()
