"""Refresh sigma and displayed G+cNG refinements after the main power run.

This runs existing diagnostics sequentially in a new directory. It neither
changes the production models nor times their performance. The main refresh
must have completed, providing the fresh 11993-node inputs and halo export.
The command prints its plan unless --execute is supplied.
"""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--cocoa-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--onecov-python", type=Path,
                        default=Path(".local/bin/python"))
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    baseline = args.baseline.resolve()
    output = args.output.resolve()
    seed = baseline/"native_seed/onecov.ini"
    commands = []

    def add(name, backend, script, *arguments):
        interpreter = getattr(args, f"{backend}_python").absolute()
        command = [str(interpreter), str(repo/"scripts"/script)]
        command.extend(str(value) for value in arguments)
        commands.append(dict(name=name, backend=backend, command=command))

    add("sigma_cocoa", "cocoa", "diagnose_sigma.py", "cocoa",
        "--inputs", baseline/"halo_onecov", "--output", output/"sigma_cocoa")
    add("sigma_onecov", "onecov", "diagnose_sigma.py", "onecov",
        "--inputs", baseline/"halo_onecov", "--cocoa", output/"sigma_cocoa",
        "--output", output/"sigma_onecov")
    add("sigma_report", "onecov", "report_sigma.py", output/"sigma_onecov",
        output/"sigma_cocoa", "--output", output/"results/sigma_abundance.json",
        "--figures", output/"results/figures")

    # Supplying Pmm also makes OneCov clamp its native hmf k support to
    # that file. Recheck finite bias normalization on this same setup.
    for nodes in (4097, 8193):
        add(f"bias_{nodes}", "onecov", "diagnose_bias.py", "onecov",
            "--config", seed, "--nodes", nodes, "--output",
            output/f"bias_onecov_{nodes}")
    add("halo_400", "onecov", "compare_halo.py", "export", seed,
        "--onecov", args.onecov.resolve(), "--mass-nodes", 400,
        "--output", output/"halo_onecov_400", "--timeout", 600)

    # The 1201-node, dz=.05 case is already in the main refresh. All new
    # SSC controls use the same 800 mass nodes, including the fine case.
    for name, radial, dz in (("300", 300, .1), ("601", 601, .1),
                             ("1201", 1201, .1),
                             ("radialfine", 2401, .05),
                             ("zfiner", 2401, .025)):
        native = output/f"ssc_projection_{name}"
        add("ssc_export_"+name, "onecov", "compare_ssc.py", "export",
            seed, "--onecov", args.onecov.resolve(), "--radial-nodes", radial,
            "--delta-z", dz, "--mass-nodes", 800, "--output", native,
            "--timeout", 600)
        add("ssc_compare_"+name, "cocoa", "compare_ssc.py", "compare",
            native, "--output", output/f"ssc_comparison_{name}")

    for name, mass, corner in (("400", 400, .001),
                               ("guard_fine", 800, .00001),
                               ("guard_finer", 800, .0000001)):
        add("trispectrum_"+name, "onecov", "compare_trispectrum.py", "export",
            seed, "--onecov", args.onecov.resolve(), "--mass-nodes", mass,
            "--corner-guard", corner, "--output", output/f"trispectrum_{name}",
            "--timeout", 600)
    add("trispectrum_guard_compare", "cocoa", "compare_trispectrum.py",
        "compare", output/"trispectrum_guard_fine", "--integration-accuracy", 2,
        "--output", output/"trispectrum_guard_cocoa")

    # Preserve the README's one-control comparisons. Their G+cNG mode ratios
    # depend on the new shared nonlinear power through the Gaussian term.
    # Each tuple is: name, k nodes, mass nodes, radial nodes, dz, corner cut.
    cases = [(f"k{n}", n, 400, 300, .5, .001) for n in (9, 17, 33, 65, 129)]
    cases.extend([
        ("radial", 65, 400, 601, .5, .001),
        ("redshift", 65, 400, 601, .25, .001),
        ("redshift_fine", 65, 400, 601, .125, .001),
        ("mass", 65, 800, 601, .25, .001),
        ("corner", 65, 400, 601, .25, .00001),
    ])
    for name, nk, nm, nr, dz, corner in cases:
        native = output/f"connected_{name}"
        add("export_"+name, "onecov", "compare_connected.py", "export",
            seed, "--onecov", args.onecov.resolve(), "--k-nodes", nk,
            "--mass-nodes", nm, "--radial-nodes", nr, "--delta-z", dz,
            "--corner-guard", corner, "--output", native, "--timeout", 600)
        add("compare_"+name, "cocoa", "compare_connected.py", "compare",
            native, "--output", output/f"connected_cocoa_{name}")
    if not args.execute:
        for command in commands:
            print(command["name"]+": "+shlex.join(command["command"]))
        return
    if output.exists():
        parser.error("choose a new output folder")
    status = json.loads((baseline/"refresh.json").read_text())
    if status["status"] != "completed":
        parser.error("the main refresh must complete before these follow-ups")
    expected = json.loads((baseline/"inputs/manifest.json").read_text())
    if expected["installed_power"]["nk"] != 11993:
        parser.error("the baseline is not the adopted 11993-node export")
    if (not 1 <= int(os.environ.get("OMP_NUM_THREADS", "0")) <= 8
            or os.environ.get("OPENBLAS_NUM_THREADS") != "1"):
        parser.error("use 1–8 OMP_NUM_THREADS and OPENBLAS_NUM_THREADS=1")
    onecov_env = os.environ.copy()
    for key in ("PYTHONPATH", "PYTHONHOME", "ROOTDIR", "VIRTUAL_ENV"):
        onecov_env.pop(key, None)
    (output/"logs").mkdir(parents=True)
    (output/"results").mkdir()
    record = dict(status="running", purpose="accuracy only",
                  baseline_sha256=sha256(baseline/"refresh.json"),
                  input_sha256=sha256(baseline/"inputs/manifest.json"),
                  scripts={p.name: sha256(p) for p in (repo/"scripts").glob("*.py")},
                  stages=commands)
    ledger = output/"refresh.json"
    try:
        for stage in commands:
            stage["status"] = "running"
            ledger.write_text(json.dumps(record, indent=2)+"\n")
            print(stage["name"], flush=True)
            env = os.environ if stage["backend"] == "cocoa" else onecov_env
            with (output/"logs"/(stage["name"]+".log")).open("w") as log:
                subprocess.run(stage["command"], cwd=repo, env=env,
                               stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=660)
            if stage["name"] == "sigma_cocoa":
                sigma = json.loads((output/"sigma_cocoa/report.json").read_text())
                if sigma["installed_power"] != expected["installed_power"]:
                    raise ValueError("sigma diagnostic changed the installed power")
            stage["status"] = "completed"
        record["status"] = "completed"
    except BaseException as error:
        record.update(status="stopped", error=str(error))
        raise
    finally:
        ledger.write_text(json.dumps(record, indent=2)+"\n")


if __name__ == "__main__":
    main()
