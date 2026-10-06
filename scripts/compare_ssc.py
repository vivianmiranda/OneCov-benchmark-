"""Compare the unbinned shear SSC projection using shared native inputs.

Export runs in OneCov's environment and calls its unchanged covELL_ssc.
Compare runs in Cocoa's environment and calls its production C kernels.
The diagnostic uses OneCov's radial nodes, matter responses, linear power
and lensing window. It tests projection, not agreement of native halo fits.
All exported distances are Mpc/h and powers/responses are (Mpc/h)^3.
"""

import argparse
import configparser
import contextlib
import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from capture_native import capture_return
from common import revision, sha256


def export(args):
    """Save a native SSC matrix and the arrays actually used to compute it."""
    sys.path.insert(0, str(args.onecov.resolve()))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_ell_space import CovELLSpace

    config = configparser.ConfigParser()
    config.read(args.source)
    if (config["observables"].getboolean("clustering")
            or config["observables"].getboolean("ggl")):
        raise ValueError("this projection comparison requires one source bin")
    config["output settings"]["directory"] = str(args.output.resolve())
    config["covariance terms"]["ssc"] = "True"
    config["covariance terms"]["nongauss"] = "False"
    config["misc"]["num_cores"] = os.environ["OMP_NUM_THREADS"]
    grid = config["covELLspace settings"]
    grid.update(ell_min="30", ell_max="3000", ell_bins="100")
    for name, value in (("integration_steps", args.radial_nodes),
                        ("delta_z", args.delta_z)):
        if value is not None:
            grid[name] = str(value)
    if args.mass_nodes is not None:
        config["halomodel evaluation"]["M_bins"] = str(args.mass_nodes)
    config_file = args.output / "onecov.ini"
    with config_file.open("w") as stream:
        config.write(stream)

    with (args.output / "native.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            started = time.perf_counter()
            terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
                Input().read_input(config_name=str(config_file)))
            tables = FileInput(bias).read_input(config_name=str(config_file))
            model = CovELLSpace(terms, obs, output, cosmo, bias, ia, hod,
                               survey, prec, tables)
            if model.n_tomo_lens != 1:
                raise ValueError("select exactly one source population")
            setup_seconds = time.perf_counter() - started
            started = time.perf_counter()
            blocks, arrays = capture_return(
                function=model.covELL_ssc,
                names=("Pmm_response", "survey_variance_mmmm", "power",
                       "ell", "sum_m_a_lm"),
                call=lambda: model.covELL_ssc(
                    bias, hod, prec, survey, obs["ELLspace"]))
            instrumented_seconds = time.perf_counter() - started

    chi = model.los_integration_chi
    # Simpson integration is a weighted sum. Integrating basis vectors
    # extracts the exact native rule; no covariance formula is reproduced.
    dchi = simpson(np.eye(len(chi)), x=chi, axis=1)
    window = model.spline_lensweight[0](chi)
    filename = args.output / "inputs.npz"
    np.savez_compressed(
        filename, chi=chi, dchi=dchi, window=window, ell=model.ellrange,
        mask_ell=arrays["ell"], mask_alm2=arrays["sum_m_a_lm"],
        power=arrays["power"], response=arrays["Pmm_response"],
        variance=arrays["survey_variance_mmmm"][0, 0],
        onecov=blocks[5][:, :, 0, 0, 0, 0, 0, 0],
        area_sr=float(survey["survey_area_lens"][0]) * (np.pi / 180)**2)
    record = dict(
        scope="unbinned shear SSC; native OneCov inputs shared with Cocoa",
        onecov=revision(args.onecov), script_sha256=sha256(__file__),
        observer_sha256=sha256(Path(__file__).with_name("capture_native.py")),
        config_sha256=sha256(config_file), inputs_sha256=sha256(filename),
        radial_nodes=len(chi), response_redshifts=len(model.los_chi),
        mass_nodes=prec["hm"]["M_bins"], shape=list(blocks[5].shape[:2]),
        setup_seconds=setup_seconds,
        instrumented_ssc_seconds=instrumented_seconds,
        timing_note="includes passive local-array capture; not a benchmark",
        omp_threads=int(os.environ["OMP_NUM_THREADS"]))
    (args.output / "export.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


def compare(args):
    """Use Cocoa's mask, shell-response and weighted-projection C routines."""
    core = args.cocoa.resolve() / "external_modules/code/cosmolike_core"
    sys.path.insert(0, str(core))
    import cosmolike_lsst_y1_interface as ci
    from cosmolike_notebook_utils.covariance.geometry import cap_mask

    metadata = json.loads((args.source / "export.json").read_text())
    filename = args.source / "inputs.npz"
    if sha256(filename) != metadata["inputs_sha256"]:
        raise ValueError("native inputs changed after export")
    data = np.load(filename)
    chi, dchi, window = data["chi"], data["dchi"], data["window"]
    area = float(data["area_sr"])
    modes = data["mask_ell"]
    if not np.array_equal(modes, np.arange(len(modes))):
        raise ValueError("the raw mask kernel requires consecutive L from zero")

    # OneCov stores sum_m |a_Lm|^2 = (2L+1) C_L. Cocoa accepts raw C_L
    # and checks its monopole against the mask's area. HEALPix's discrete
    # cap differs slightly from the nominal survey area. Use its monopole
    # area in the kernel, then restore OneCov's nominal-area normalization.
    mask = np.ascontiguousarray(data["mask_alm2"] / (2*modes+1))
    pixel_area = np.sqrt(4*np.pi*mask[0])
    backend = ci.covariance
    variance = backend.covariance_ssc_mask_variance(
        mask_cl=mask, area_sr=pixel_area, distance=chi, power=data["power"])
    variance *= (pixel_area/area)**2
    pair = np.ascontiguousarray(np.broadcast_to(window**2,
                                               data["response"].T.shape))
    shell = backend.covariance_ssc_shell_response(
        distance=chi, signal=np.zeros(len(data["ell"])), pair_window=pair,
        mean_window=np.zeros_like(pair),
        power_response=np.ascontiguousarray(data["response"].T))
    cocoa = backend.covariance_project(
        left=shell, right=shell, weight=np.ascontiguousarray(dchi*variance))

    # Separately replace only the pixelized footprint by Cocoa's analytic
    # spherical cap, with the same nominal area and multipole cutoff.
    analytic_variance = backend.covariance_ssc_mask_variance(
        mask_cl=cap_mask(area, len(modes)-1), area_sr=area,
        distance=chi, power=data["power"])
    analytic = backend.covariance_project(
        left=shell, right=shell,
        weight=np.ascontiguousarray(dchi*analytic_variance))
    onecov = data["onecov"]
    rms = np.sqrt(np.diag(onecov))
    residual = (cocoa-onecov)/rms[:, None]/rms[None, :]
    mask_change = (analytic-onecov)/rms[:, None]/rms[None, :]
    variance_error = np.max(np.abs(variance*chi**2/data["variance"]-1))
    eigenvalues = np.linalg.eigvalsh(onecov/rms[:, None]/rms[None, :])
    record = dict(
        scope=metadata["scope"], export=metadata,
        max_variance_scaled_difference=float(np.max(np.abs(residual))),
        mask_variance_max_fractional_difference=float(variance_error),
        analytic_cap_max_variance_scaled_change=float(np.max(np.abs(mask_change))),
        pixel_area_fractional_difference=float(pixel_area/area-1),
        minimum_simpson_weight=float(dchi.min()),
        minimum_ssc_correlation_eigenvalue=float(eigenvalues.min()),
        max_ssc_correlation_asymmetry=float(np.max(
            np.abs(onecov-onecov.T)/rms[:, None]/rms[None, :])),
        interface_sha256=sha256(ci.__file__), script_sha256=sha256(__file__),
        passed=bool(np.max(np.abs(residual)) < 1.e-11
                    and eigenvalues.min() > -1.e-10),
        remaining="independent halo responses, radial and interpolation convergence")
    np.savez_compressed(args.output / "matrices.npz", onecov=onecov,
                        cocoa=cocoa, analytic_cap=analytic, ell=data["ell"])
    (args.output / "comparison.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))
    if not record["passed"]:
        raise SystemExit("SSC projection differs; results preserved")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("export", "compare"))
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--cocoa", type=Path,
                        default=Path(__file__).resolve().parents[2]/"cocoa/Cocoa")
    parser.add_argument("--radial-nodes", type=int)
    parser.add_argument("--delta-z", type=float)
    parser.add_argument("--mass-nodes", type=int)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output directory")
    if int(os.environ.get("OMP_NUM_THREADS", "0")) < 1:
        parser.error("set OMP_NUM_THREADS explicitly")
    args.output.mkdir(parents=True)
    signal.alarm(args.timeout)
    if args.mode == "export":
        export(args)
    else:
        compare(args)


if __name__ == "__main__":
    main()
