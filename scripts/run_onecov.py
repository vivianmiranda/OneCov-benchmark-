"""Run a bounded native OneCovariance LSST Y1 Fourier-space pilot.

The generated INI describes one source bin, optionally with two galaxy
lens bins. G, G+SSC and G+cNG are separate cases. Shared C_ell input isolates
Gaussian assembly; native spectra and shared P(k,z) test different parts
of the pipeline and must be reported separately.

The script validates exported input hashes, writes the native INI and
uses OneCovariance's own parser and CLI. It never patches numerical code.
Each run gets a new directory so a failed or timed-out run is inspectable.
Set OMP_NUM_THREADS before starting; the runner fixes BLAS to one thread.
"""

import argparse
import configparser
import contextlib
import importlib.metadata
import json
import os
import platform
import resource
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from common import revision, sha256


def make_config(template, inputs, output, manifest, args, workers):
    """Combine numerical settings with the exported LSST physical inputs.

    Arguments:
        template = INI containing pilot numerical choices.
        inputs, output = absolute input and result directories.
        manifest = metadata from prepare_lsst_y1.py.
        args = parsed case, component and spectra-source selections.
        workers = worker count read from OMP_NUM_THREADS.
    Returns:
        ConfigParser for the native OneCovariance CLI; no files are written.
    """
    config = configparser.ConfigParser()
    config.read(filenames=template)
    # Components stay separate even when a copied template has different
    # output switches. The --terms option adds one contribution to G.
    config["covariance terms"]["gauss"] = "True"
    config["covariance terms"]["split_gauss"] = "True"
    config["covariance terms"]["ssc"] = str(args.terms == "ssc")
    config["covariance terms"]["nongauss"] = str(args.terms == "connected")
    config["observables"]["ggl"] = str(args.case == "3x2")
    config["observables"]["clustering"] = str(args.case == "3x2")
    config["bias"]["bias_files"] = str(inputs / "bias.txt")
    config["misc"] = {"num_cores": str(workers)}

    # Narrow bands keep the integer-grid assembly test small. They change
    # the measured estimator, not the covariance formula or the spectra.
    if args.band_limits is not None:
        lower, upper = args.band_limits
        for tracer in ("lensing", "clustering"):
            config["covELLspace settings"][f"ell_min_{tracer}"] = str(lower)
            config["covELLspace settings"][f"ell_max_{tracer}"] = str(upper)

    cosmo = manifest["cosmology"]
    density_strings = []
    for value in manifest["lens_density_arcmin2"]:
        density_strings.append(str(value))
    config["cosmo"] = {
        "sigma8": str(manifest["sigma8"]),
        "h": str(cosmo["H0"] / 100.0),
        "omega_m": str(cosmo["omegam"]),
        "omega_b": str(cosmo["omegab"]),
        "omega_de": str(1.0 - cosmo["omegam"]),
        "w0": str(cosmo["w"]),
        "wa": str(cosmo["w0pwa"] - cosmo["w"]),
        "ns": str(cosmo["ns"]),
        "m_nu": "0",
        "neff": "3.044",
    }
    config["survey specs"] = {
        "survey_area_lensing_in_deg2": str(manifest["area_deg2"]),
        "survey_area_clust_in_deg2": str(manifest["area_deg2"]),
        "survey_area_ggl_in_deg2": str(manifest["area_deg2"]),
        "n_eff_lensing": str(manifest["source_density_arcmin2"]),
        "n_eff_clust": ", ".join(density_strings),
        "ellipticity_dispersion": str(manifest["sigma_e_component"]),
    }
    config["redshift"] = {
        "z_directory": str(inputs),
        "zlens_file": "source.nz",
        "zclust_file": "lens.nz",
        "value_loc_in_lensbin": "mid",
        "value_loc_in_clustbin": "mid",
    }

    # Each diagnostic supplies only the ingredient named by --spectra.
    # A shared Pmm does not replace OneCov's halo-model linear power or
    # responses. Supplying Cells bypasses their projection, not halo setup.
    config["tabulated inputs files"] = {}
    if args.spectra == "shared-cells":
        config["tabulated inputs files"] = {
            "Cell_directory": str(inputs),
            "Cmm_file": "Cmm.txt",
        }
        if args.case == "3x2":
            config["tabulated inputs files"]["Cgg_file"] = "Cgg.txt"
            config["tabulated inputs files"]["Cgm_file"] = "Cgm.txt"
    elif args.spectra == "shared-power":
        config["tabulated inputs files"] = {
            "powspec_directory": str(inputs),
            "Pmm_file": "Pmm.txt",
        }
    config["output settings"]["directory"] = str(output)
    return config


def check_native_inputs(config_file, onecov, args):
    """Ask OneCov's parser to read the actual files before numerical work.

    Arguments:
        config_file = generated INI path.
        onecov = source checkout, added to the Python import path.
        args = requested case and supplied-spectrum choice.
    Returns:
        Actual CAMB version/path, useful when imported from a source checkout.
        Parser messages go to the caller's redirected log.
    Raises:
        RuntimeError if a requested table was ignored instead of loaded.
        OneCov's own input exceptions are preserved for diagnosis.
    """
    sys.path.insert(0, str(onecov))
    import camb
    from onecov.cov_input import FileInput, Input

    parsed = Input().read_input(config_name=str(config_file))
    bias = parsed[4]
    tables = FileInput(bias).read_input(config_name=str(config_file))
    if args.case == "3x2" and tables["zet_dep_bias"] is None:
        raise RuntimeError("OneCov ignored the bias table; do not fall back to HOD")
    if args.spectra == "shared-cells":
        if tables["Cxy"]["mm"] is None:
            raise RuntimeError("OneCov ignored Cmm; check the explicit four-column file")
        if args.case == "3x2":
            if tables["Cxy"]["gg"] is None or tables["Cxy"]["gm"] is None:
                raise RuntimeError("OneCov ignored a galaxy C_ell table")
    if args.spectra == "shared-power" and tables["Pxy"]["mm"] is None:
        raise RuntimeError("OneCov ignored Pmm; check its z-major table format")
    return {"version": camb.__version__, "path": camb.__file__}


def main():
    """Write a reproducible run directory and execute one time-limited pilot.

    Returns:
        None. Writes native inputs, logs and run.json with completion status.
    Raises:
        SystemExit on configuration, dependency, timeout or native-run errors.
        Partial output is retained; select a new output directory to retry.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=("shear", "3x2"), default="shear")
    parser.add_argument("--terms", choices=("gaussian", "ssc", "connected"),
                        default="gaussian")
    parser.add_argument("--spectra", choices=("shared-cells", "shared-power", "native"),
                        default="shared-cells")
    parser.add_argument("--template", type=Path,
                        default=Path(__file__).resolve().parents[1] /
                        "configs" / "onecov_pilot.ini")
    parser.add_argument("--timeout", type=float, default=600.0,
                        help="wall-time budget in seconds, default 600")
    parser.add_argument("--band-limits", type=int, nargs=2,
                        metavar=("LOWER", "UPPER"),
                        help="override output band range [LOWER, UPPER)")
    parser.add_argument("--prepare-only", action="store_true",
                        help="write INI and provenance without importing OneCov")
    args = parser.parse_args()
    onecov = args.onecov.resolve()
    inputs = args.inputs.resolve()
    output = args.output.resolve()
    if not (onecov / "covariance.py").is_file():
        parser.error(f"{onecov}: missing native covariance.py; check --onecov")
    if not args.template.is_file():
        parser.error(f"{args.template}: missing numerical INI template")
    if output.exists():
        parser.error(f"{output} exists; select a new output directory")
    if not 0 < args.timeout < float("inf"):
        parser.error("--timeout must be a finite positive number of seconds")
    if args.band_limits is not None:
        lower, upper = args.band_limits
        if lower < 2 or upper <= lower:
            parser.error("--band-limits requires 2 <= LOWER < UPPER")
    try:
        workers = int(os.environ.get("OMP_NUM_THREADS", "0"))
    except ValueError:
        parser.error("set OMP_NUM_THREADS to a single positive integer")
    if workers < 1:
        parser.error("set OMP_NUM_THREADS explicitly, e.g. export OMP_NUM_THREADS=8")

    # Set these before importing any numerical packages in the preflight
    # or child. Only OneCov's chosen parallel work may use multiple cores.
    for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[variable] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    manifest_file = inputs / "manifest.json"
    manifest = json.loads(manifest_file.read_text())
    for name, expected in manifest["files"].items():
        if sha256(filename=inputs / name) != expected:
            parser.error(f"{inputs / name} changed; create a fresh input export")

    config = make_config(template=args.template, inputs=inputs, output=output,
                         manifest=manifest, args=args, workers=workers)
    output.mkdir(parents=True)
    os.environ["MPLCONFIGDIR"] = str(output / "matplotlib")
    config_file = output / "onecov.ini"
    with config_file.open(mode="w") as stream:
        config.write(fp=stream)
    (output / "input_manifest.json").write_text(manifest_file.read_text())

    command = [sys.executable, str(onecov / "covariance.py"), str(config_file)]
    report = {
        "status": "prepared",
        "case": args.case,
        "terms": args.terms,
        "spectra": args.spectra,
        "onecov": revision(directory=onecov),
        "benchmark": revision(directory=Path(__file__).resolve().parents[1]),
        "input_manifest_sha256": sha256(filename=manifest_file),
        "config_sha256": sha256(filename=config_file),
        "runner_sha256": sha256(filename=__file__),
        "template_sha256": sha256(filename=args.template),
        "created_utc": datetime.now(tz=timezone.utc).isoformat(),
        "command": command,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "workers": workers,
        "blas_threads": 1,
        "timeout_seconds": args.timeout,
        "packages": {},
    }
    for name in ("numpy", "scipy", "astropy", "hmf", "healpy", "camb", "levin"):
        try:
            report["packages"][name] = importlib.metadata.version(distribution_name=name)
        except importlib.metadata.PackageNotFoundError:
            report["packages"][name] = None
    report_file = output / "run.json"
    report_file.write_text(json.dumps(report, indent=2) + "\n")
    if args.prepare_only:
        print(f"Prepared {config_file}; no OneCov numerical work was run.")
        return

    # Parsing has no covariance work. Record failures separately so a failed
    # import or an ignored table can never become an apparent fast timing.
    try:
        with (output / "preflight.log").open(mode="w") as log:
            with contextlib.redirect_stdout(log):
                report["camb_runtime"] = check_native_inputs(
                    config_file=config_file, onecov=onecov, args=args,
                )
    except Exception as error:
        report["status"] = "preflight_failed"
        report["error"] = str(error)
        report_file.write_text(json.dumps(report, indent=2) + "\n")
        raise

    print(f"Running {args.case}, {args.terms}, {args.spectra}", flush=True)
    print(f"Native log: {output / 'native.log'}", flush=True)
    started = time.perf_counter()
    with (output / "native.log").open(mode="w") as log:
        process = subprocess.Popen(args=command, cwd=output, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            returncode = process.wait(timeout=args.timeout)
            report["status"] = "completed" if returncode == 0 else "failed"
        except subprocess.TimeoutExpired:
            # Stop the whole child group: a time budget must not leave any
            # worker processes consuming cores after the parent exits.
            os.killpg(process.pid, signal.SIGTERM)
            try:
                returncode = process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                returncode = process.wait()
            report["status"] = "timed_out"
    report["returncode"] = returncode
    report["native_cli_wall_seconds"] = time.perf_counter() - started
    report["timing_scope"] = "fresh native CLI: imports, setup, computation and output"
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    rss_scale = 1024.0
    if sys.platform == "darwin":
        rss_scale = 1.0
    report["peak_child_rss_mib"] = usage.ru_maxrss * rss_scale / 1024**2
    report["memory_scope"] = "largest child RSS; not a sum across simultaneous workers"
    report_file.write_text(json.dumps(report, indent=2) + "\n")
    print(f"{report['status']}: {report['native_cli_wall_seconds']:.2f} s")
    if report["status"] != "completed":
        raise SystemExit(f"Inspect {output / 'native.log'}; partial files retained")


if __name__ == "__main__":
    main()
