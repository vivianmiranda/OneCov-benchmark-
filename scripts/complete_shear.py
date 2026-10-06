"""Compute a complete one-source covariance in each code, separately.

The observable is the LSST Y1 source-bin-3 shear power at ell=30,60,...,
3000. Both calculations use the band-centre approximation with width 30:
G has 30 times the single-multipole mode count; SSC and cNG are evaluated
at the band centres. This avoids the codes' different broad-band weights.
With --space real, instead compute xi+ and xi- in eight logarithmic bins
from 2.5 to 250 arcminutes, including their cross-covariance. CoCoA uses
full-sky spherical annuli; OneCov uses flat-sky Bessel annuli and the
numerical settings in its shipped config_3x2pt_rcf.ini. Neither pilot is
the full LSST covariance. Refinements must be requested explicitly.

The Gaussian projections receive the same CAMB nonlinear power. Each code
keeps its own halo model, linear-power preparation, response and footprint.
No numerical source or callable is replaced. Run the two environments
sequentially and choose a new output folder for every calculation.
"""

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np

from common import revision, sha256
from run_onecov import make_config


def theta_matrix(blocks):
    """Join native ++, +-, -- blocks, ordered xi+ bins then xi- bins.

    Each block has two angular axes and six singleton population/tomography
    axes. The transpose fills the reversed observable ordering (-,+), not
    an average of independently computed entries within either auto block.
    """
    plus, cross, minus = [block[:, :, 0, 0, 0, 0, 0, 0] for block in blocks]
    return np.block([[plus, cross], [cross.T, minus]])


def onecov(args, manifest):
    """Call native G/SSC/cNG for Fourier band centres or real-space annuli.

    args supplies input/output paths and numerical refinement controls.
    manifest supplies the exported cosmology, density and source bin.
    Return dimensionless covariance components and run provenance. Fourier
    matrices are [100,100]; real-space matrices are [16,16], xi+ then xi-.
    """
    sys.path.insert(0, str(args.onecov.resolve()))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_ell_space import CovELLSpace
    from onecov.cov_theta_space import CovTHETASpace

    # The existing input adapter sets the survey and shared P(k,z) file.
    # Both non-Gaussian switches must be enabled with split_gauss=True,
    # otherwise OneCov can combine SSC with its connected output.
    choices = argparse.Namespace(
        terms="connected", case="shear", spectra="shared-power",
        band_limits=None,
    )
    template = Path(__file__).resolve().parents[1]/"configs/onecov_pilot.ini"
    if args.space == "real":
        template = args.onecov.resolve()/"config_files/config_3x2pt_rcf.ini"
    config = make_config(
        template=template,
        inputs=args.inputs.resolve(), output=args.output.resolve(),
        manifest=manifest, args=choices,
        workers=int(os.environ["OMP_NUM_THREADS"]),
    )
    config["covariance terms"]["ssc"] = "True"
    grid = config["covELLspace settings"]
    if args.space == "fourier":
        grid.update(ell_min="30", ell_max="3000", ell_bins="100", ell_type="lin")
    # Omit broad output bands. OneCov's public unbinned G routine uses
    # delta_ell=30 from this linear grid; its SSC/cNG retain both ell axes.
    for tracer in ("lensing", "clustering"):
        for key in ("min", "max", "bins", "type"):
            config.remove_option("covELLspace settings", f"ell_{key}_{tracer}")
    if args.space == "real":
        config["observables"]["est_shear"] = "xi_pm"
        config["IA"]["A_IA"] = "0"
        grid["mult_shear_bias"] = "0"
        # Copy the native real-space example's numerical controls. Only
        # the observed angular bins and physical IA/calibration inputs
        # change; do not tune OneCov's integrator to CoCoA's cutoff.
        config["covTHETAspace settings"].update(
            theta_min="2.5", theta_max="250", theta_bins="8",
            theta_type="log", theta_min_lensing="2.5",
            theta_max_lensing="250", theta_bins_lensing="8",
            theta_type_lensing="log", xi_pp="True", xi_mm="True",
        )
    overrides = (
        ("covELLspace settings", "integration_steps", args.radial_nodes),
        ("covELLspace settings", "delta_z", args.delta_z),
        ("covELLspace settings", "tri_delta_z", args.tri_delta_z),
        ("halomodel evaluation", "M_bins", args.mass_nodes),
        ("trispec evaluation", "log10k_bins", args.k_nodes),
        ("covELLspace settings", "ell_bins", args.ell_nodes),
        ("covTHETAspace settings", "theta_accuracy", args.theta_accuracy),
    )
    for section, key, value in overrides:
        if value is not None:
            config[section][key] = str(value)
    config["output settings"]["corrmatrix_plot"] = "False"
    filename = args.output/"onecov.ini"
    with filename.open("w") as stream:
        config.write(stream)

    started = time.perf_counter()
    terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
        Input().read_input(config_name=str(filename)))
    tables = FileInput(bias).read_input(config_name=str(filename))
    if tables["Pxy"]["mm"] is None:
        raise ValueError("OneCov did not load the supplied nonlinear power")
    constructor = CovTHETASpace if args.space == "real" else CovELLSpace
    model = constructor(terms, obs, output, cosmo, bias, ia, hod,
                        survey, prec, tables)
    ell = np.arange(30, 3001, 30, dtype=float)
    if args.space == "fourier":
        np.testing.assert_array_equal(model.ellrange, ell)
    else:
        np.testing.assert_allclose(model.theta_ul_bins_lensing,
                                   np.geomspace(2.5, 250, 9), rtol=1.e-13)
    if model.n_tomo_lens != 1 or model.gg or model.gm:
        raise ValueError("the complete pilot requires one shear source bin")
    if model.ellrange_lensing_ul is not None:
        raise ValueError("the pilot must use the stated band-centre estimator")
    setup_seconds = time.perf_counter()-started

    started = time.perf_counter()
    calculate = model.calc_covTHETA if args.space == "real" else model.calc_covELL
    gaussian, connected, ssc = calculate(
        obs, output, bias, hod, survey, prec, tables,
    )
    construction_seconds = time.perf_counter()-started
    # Native output has six tracer groups; mmmm is group 5. Split G
    # has three entries per group: signal, mixed signal/noise, pure noise.
    # Singleton population and tomography axes select this source only.
    matrices = {}
    if args.space == "fourier":
        matrices = {
            "gaussian": sum(gaussian[15:18])[:, :, 0, 0, 0, 0, 0, 0],
            "ssc": ssc[5][:, :, 0, 0, 0, 0, 0, 0],
            "cng": connected[5][:, :, 0, 0, 0, 0, 0, 0],
            "ell": ell,
            "signal": model.Cell_mm[:, 0, 0, 0],
        }
    if args.space == "real":
        # Real-space groups 7,8,9 are ++,+-,--. The split Gaussian list
        # has signal, mixed and pure-noise pieces for each of these groups.
        matrices = {
            "gaussian": theta_matrix([sum(gaussian[21:24]),
                                       sum(gaussian[24:27]),
                                       sum(gaussian[27:30])]),
            "ssc": theta_matrix(ssc[7:10]),
            "cng": theta_matrix(connected[7:10]),
            "theta_edges_arcmin": model.theta_ul_bins_lensing,
            "signal": np.concatenate([model.xi_plus.ravel(),
                                      model.xi_minus.ravel()]),
        }
    record = {
        "onecov": revision(args.onecov),
        "config_sha256": sha256(filename),
        "template": str(template.relative_to(args.onecov.resolve()))
                    if args.space == "real" else str(template),
        "template_sha256": sha256(template),
        "setup_seconds": setup_seconds,
        "construction_seconds": construction_seconds,
        "radial_nodes": len(model.los_integration_chi),
        "k_nodes": config.getint("trispec evaluation", "log10k_bins"),
        "mass_nodes": config.getint("halomodel evaluation", "M_bins"),
        "delta_z": grid.getfloat("delta_z"),
        "tri_delta_z": grid.getfloat("tri_delta_z"),
        "mass_min": 10**config.getfloat("halomodel evaluation", "log10M_min"),
        "ell_samples": model.ellrange.tolist(),
        "theta_accuracy": obs["THETAspace"]["theta_acc"],
    }
    if args.space == "real":
        record["transform_ell_endpoints"] = model.ell_fourier_integral[[0,-1]].tolist()
        record["integration_intervals"] = model.integration_intervals
    return matrices, record


def cocoa(args, manifest):
    """Call the production forecast for the same source and observables.

    The production call uses one-multipole operators. Dividing only G by
    30 matches OneCov's explicit band-centre mode count. SSC/cNG correlate
    the band centres through shared matter modes and receive no such factor.
    The numerical settings, native predictions and timings remain recorded.
    Real-space instead requests both spin kernels and keeps G unchanged:
    its annular bin averaging and noise are already included in the result.
    """
    runtime = args.cocoa.resolve()
    project = runtime/"projects/lsst_y1"
    core = runtime/"external_modules/code/cosmolike_core"
    sys.path.insert(0, str(project/"covariance"))
    sys.path.insert(0, str(project/"interface"))
    sys.path.insert(0, str(core))
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize, compute
    from cosmolike_notebook_utils.covariance.forecast import _json_array

    settings = configuration(
        accuracy_boost=args.accuracy_boost,
        integration_accuracy=args.integration_accuracy,
        gaussian={"nonlimber": False, "ia": "none"},
    )
    if settings["cosmology"] != manifest["cosmology"]:
        raise ValueError("current LSST cosmology differs from the input export")
    source = manifest["source_bin_1based"]-1
    expected = (
        (settings["area_deg2"], manifest["area_deg2"]),
        (settings["source_density_arcmin2"][source],
         manifest["source_density_arcmin2"]),
        (settings["sigma_e_component"][source], manifest["sigma_e_component"]),
    )
    for actual, reference in expected:
        if actual != reference:
            raise ValueError("survey density, area or shape noise differs")
    if sha256(project/settings["source_file"]) != manifest["source_nz_sha256"]:
        raise ValueError("LSST source distribution differs from the export")
    ell = np.arange(30, 3001, 30, dtype=np.int32)
    settings["band_first"] = ell
    settings["band_last"] = ell.copy()
    if args.space == "real":
        settings["theta_edges_arcmin"] = np.geomspace(2.5, 250, 9)
    started = time.perf_counter()
    tables = initialize(interface=ci, settings=settings)
    setup_seconds = time.perf_counter()-started

    # Verify the shared power against the actual CAMB arrays installed in
    # CoCoA. OneCov receives the z-major ASCII subset of these same arrays.
    supplied = np.loadtxt(args.inputs/"Pmm.txt")
    wave = 10.0**tables["log10k_2D"]
    redshift = tables["z_2D"]
    power = np.exp(tables["lnP_nonlinear"].reshape(
        (len(redshift), len(wave)), order="F"))
    nredshift = len(supplied)//len(wave)
    np.testing.assert_allclose(supplied[:, 0],
                               np.repeat(redshift[:nredshift], len(wave)),
                               rtol=1.e-13, atol=0)
    np.testing.assert_allclose(supplied[:, 1], np.tile(wave, nredshift),
                               rtol=1.e-13, atol=0)
    np.testing.assert_allclose(supplied[:, 2], power[:nredshift].ravel(),
                               rtol=1.e-11, atol=0)

    # The initialized core retains all lens/source bins for crossed spectra.
    # Select only source bin 3: one E-mode row in Fourier space, or the
    # two observed correlations xi+ and xi- in real space.
    field = len(settings["lens_density_arcmin2"])+source
    rows = np.array([[0, field, field]], dtype=np.int32)
    if args.space == "real":
        rows = np.array([[0, field, field], [1, field, field]], dtype=np.int32)

    def progress(stage, seconds):
        print(f"{stage}: {seconds:.2f} s", flush=True)

    result = compute(
        interface=ci, settings=settings, space=args.space, rows=rows,
        backend=ci.covariance, progress=progress,
    )
    matrices = {
        "gaussian": result["gaussian"]/30.0,
        "ssc": result["ssc"],
        "cng": result["cng"],
        "ell": ell.astype(float),
        "signal": result["signal"].ravel(),
    }
    if args.space == "real":
        matrices["gaussian"] = result["gaussian"]
        del matrices["ell"]
        matrices["theta_edges_arcmin"] = settings["theta_edges_arcmin"]
    record = {
        "core": revision(core),
        "project": revision(project),
        "interface_sha256": sha256(ci.__file__),
        "setup_seconds": setup_seconds,
        "construction_seconds": result["stages_s"]["total"],
        "stages_seconds": result["stages_s"],
        "settings": json.loads(json.dumps(settings, default=_json_array)),
        "mass_min": float(np.exp(settings["lnm_edges"][0])),
    }
    return matrices, record


def main():
    """Run one backend and archive all components without repairing modes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("backend", choices=("onecov", "cocoa"))
    parser.add_argument("--space", choices=("fourier", "real"), default="fourier")
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--cocoa", type=Path, default=Path("../cocoa/Cocoa"))
    parser.add_argument("--k-nodes", type=int)
    parser.add_argument("--mass-nodes", type=int)
    parser.add_argument("--radial-nodes", type=int)
    parser.add_argument("--delta-z", type=float)
    parser.add_argument("--tri-delta-z", type=float)
    parser.add_argument("--ell-nodes", type=int)
    parser.add_argument("--theta-accuracy", type=float)
    parser.add_argument("--accuracy-boost", type=int, default=1)
    parser.add_argument("--integration-accuracy", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--accuracy-only", action="store_true",
                        help="omit timings when other numerical jobs are running")
    args = parser.parse_args()
    if args.space == "fourier":
        defaults = dict(k_nodes=129, mass_nodes=400, radial_nodes=601,
                        delta_z=0.05, tri_delta_z=0.125)
        for key, value in defaults.items():
            if getattr(args, key) is None:
                setattr(args, key, value)
    if args.output.exists():
        parser.error("choose a new output directory")
    if int(os.environ.get("OMP_NUM_THREADS", "0")) < 1:
        parser.error("set OMP_NUM_THREADS explicitly")
    manifest_path = args.inputs/"manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest["source_bin_1based"] != 3:
        parser.error("this pilot uses original LSST source bin 3")
    for filename, fingerprint in manifest["files"].items():
        if sha256(args.inputs/filename) != fingerprint:
            raise ValueError(f"exported input changed: {filename}")
    signal.alarm(args.timeout)
    args.output.mkdir(parents=True)
    if args.backend == "onecov":
        matrices, record = onecov(args=args, manifest=manifest)
    else:
        matrices, record = cocoa(args=args, manifest=manifest)

    # Concurrent jobs are suitable for numerical comparisons, but their
    # elapsed times cannot establish either code's performance. Keep those
    # times out of the result consumed by the cross-code report and plots.
    if args.accuracy_only:
        for key in ("setup_seconds", "construction_seconds", "stages_seconds"):
            record.pop(key, None)
        record["timing_note"] = "Accuracy-only run; concurrent validation job."

    matrices["total"] = matrices["gaussian"]+matrices["ssc"]+matrices["cng"]
    # Preserve the native output even if a subsequent diagnostic fails.
    # No matrix is symmetrized, masked or eigenvalue-repaired in this archive.
    filename = args.output/"covariance.npz"
    np.savez_compressed(filename, **matrices)
    ndata = 16 if args.space == "real" else 100
    rms = np.sqrt(np.diag(matrices["total"]))
    asymmetries = {}
    for name in ("gaussian", "ssc", "cng", "total"):
        matrix = matrices[name]
        if matrix.shape != (ndata, ndata) or not np.isfinite(matrix).all():
            raise ValueError(f"invalid {name} covariance; inspect the native log")
        asymmetry = np.max(np.abs(matrix-matrix.T)/rms[:, None]/rms[None, :])
        asymmetries[name] = float(asymmetry)
        if asymmetry > 1.e-4:
            raise ValueError(f"asymmetric {name} covariance: {asymmetry}")
    normalized = matrices["total"]/rms[:, None]/rms[None, :]
    # A double numerical projection can leave small antisymmetric residuals.
    # Eigenvalues describe its symmetric part; the raw residual is reported.
    minimum = float(np.linalg.eigvalsh((normalized+normalized.T)/2)[0])
    if not np.isfinite(minimum) or minimum <= 0:
        raise ValueError("total covariance is not positive definite")
    record.update(
        backend=args.backend,
        scope="one LSST Y1 source bin; Fourier band-centre approximation",
        source_bin_1based=3, delta_ell=30,
        inputs=manifest, input_manifest_sha256=sha256(manifest_path),
        script_sha256=sha256(__file__), covariance_sha256=sha256(filename),
        omp_threads=int(os.environ["OMP_NUM_THREADS"]),
        minimum_total_correlation_eigenvalue=minimum,
        space=args.space, normalized_asymmetry=asymmetries,
    )
    if args.space == "real":
        record["scope"] = "one LSST Y1 source bin; xi+, xi- and their cross"
        record["delta_ell"] = None
        record["geometry"] = "full sky" if args.backend == "cocoa" else "flat sky"
    (args.output/"report.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"Saved all four {ndata}x{ndata} matrices to {args.output}", flush=True)


if __name__ == "__main__":
    main()
