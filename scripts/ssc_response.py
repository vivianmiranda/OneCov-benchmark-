"""Start the SSC comparison with shared halo-response ingredients.

Export runs in the OneCov environment; compare runs in the Cocoa environment.
No angular SSC matrix is computed. The test isolates the response formula
from the halo mass function, concentration, mask and radial projection.
All powers and dimensional responses use (Mpc/h)^3; k uses h/Mpc.
"""

import argparse
import configparser
import contextlib
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline

from common import revision, sha256


def export_response(source, output, onecov):
    """Save OneCov's native response and its ingredients at z=0, 0.5 and 1.

    Arguments:
        source = prepared one-source OneCov INI, with the LSST cosmology.
        output = new directory; onecov = unmodified source checkout.
    Returns:
        None. Writes arrays, numerical settings, provenance and a native log.
    """
    sys.path.insert(0, str(onecov))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_polyspectra import PolySpectra

    config = configparser.ConfigParser()
    config.read(filenames=source)
    if (config["observables"].getboolean("clustering")
            or config["observables"].getboolean("ggl")):
        raise ValueError("this first response test requires shear only")
    config["output settings"]["directory"] = str(output)
    config["covariance terms"]["ssc"] = "True"
    config_file = output / "onecov.ini"
    with config_file.open(mode="w") as stream:
        config.write(stream)

    arrays = {"redshift": np.array([0.0, 0.5, 1.0])}
    with (output / "native.log").open(mode="w") as log:
        with contextlib.redirect_stdout(log):
            terms, obs, settings, cosmo, bias, ia, hod, survey, prec = (
                Input().read_input(config_name=str(config_file)))
            tables = FileInput(bias).read_input(config_name=str(config_file))
            model = PolySpectra(0.0, terms, obs, cosmo, bias, hod,
                                survey, prec, tables)
            for row, redshift in enumerate(arrays["redshift"]):
                model.update_mass_func(zet=redshift, bias_dict=bias,
                                       hod_dict=hod, prec=prec)
                response = model.powspec_responses(
                    bias_dict=bias, hod_dict=hod, hm_prec=prec["hm"])[2][:, 0]
                i11 = model.halo_model_integral_I_alpha_x(
                    bias_dict=bias, hod_dict=hod, hm_prec=prec["hm"],
                    alpha=1, type_x="m")[:, 0]
                full_i12 = model.halo_model_integral_I_alpha_xy(
                    bias_dict=bias, hod_dict=hod, hm_prec=prec["hm"],
                    alpha=1, type_x="m", type_y="m")
                i12 = np.diagonal(full_i12[:, :, 0, 0])

                # OneCov's I_alpha_xy handles biased moments only. Its
                # existing one-halo helper provides the unbiased I02.
                # Calling it exposes an ingredient; no method is replaced.
                i02 = model._PolySpectra__P_xy_1h(
                    bias_dict=bias, hod_dict=hod, hm_prec=prec["hm"],
                    type_x="m", type_y="m")[:, 0]
                k = model.mass_func.k.copy()
                linear = model.mass_func.power.copy()

                # Reproduce the slope used INSIDE powspec_responses.
                # dln(k^3 P)/dlnk = 3 + dlnP/dlnk. Cocoa's low-level
                # routine expects the latter, so subtract three explicitly.
                smoothing = 5 * np.min(np.diff(k))
                derivative = UnivariateSpline(
                    x=np.log(k), y=np.log(k**3 * linear), k=3,
                    s=smoothing).derivative()
                slope = derivative(np.log(k)) - 3.0
                # A separate diagnostic adds the derivative of I11^2,
                # as in the corrected Takada-Hu two-halo prescription.
                i11_slope = UnivariateSpline(
                    x=np.log(k), y=np.log(i11), k=3, s=0).derivative()(np.log(k))

                if row == 0:
                    arrays["k_h_Mpc"] = k
                    for name in ("linear", "i11", "i02", "i12", "slope",
                                 "two_halo_slope", "response_onecov"):
                        arrays[name] = np.empty(shape=(3, len(k)))
                for name, value in (("linear", linear), ("i11", i11),
                                    ("i02", i02), ("i12", i12),
                                    ("slope", slope),
                                    ("two_halo_slope", slope + 2 * i11_slope),
                                    ("response_onecov", response)):
                    arrays[name][row] = value
    filename = output / "ingredients.npz"
    np.savez_compressed(file=filename, **arrays)
    record = {
        "scope": "native OneCov matter response and shared halo ingredients",
        "onecov": revision(directory=onecov),
        "source_config_sha256": sha256(filename=source),
        "resolved_config_sha256": sha256(filename=config_file),
        "ingredients_sha256": sha256(filename=filename),
        "script_sha256": sha256(filename=__file__),
        "omp_threads": int(os.environ["OMP_NUM_THREADS"]),
        "power_and_response_units": "(Mpc/h)^3",
        "moment_note": "Native matter moments temporarily extend log10M_min to 2; configured mass limits alone do not describe their integration range",
    }
    (output / "export.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Exported native SSC responses at three redshifts to {output}")


def compare_response(source, output):
    """Test Cocoa's existing response kernel with OneCov's supplied moments.

    The first result matches OneCov's slope and disables fractional
    transfer to another power spectrum. A second result changes only
    the differentiated spectrum to I11^2 P_linear. This second result is
    a prescription diagnostic, not a native Cocoa SSC prediction.
    """
    import cosmolike_lsst_y1_interface as ci

    metadata = json.loads((source / "export.json").read_text())
    filename = source / "ingredients.npz"
    if sha256(filename=filename) != metadata["ingredients_sha256"]:
        raise ValueError("SSC ingredients changed after export")
    data = np.load(file=filename)
    shape = data["linear"].shape
    halo_power = data["i11"]**2 * data["linear"] + data["i02"]
    inputs = np.ascontiguousarray([
        data["linear"].ravel(), halo_power.ravel(), data["i11"].ravel(),
        data["i02"].ravel(), data["i12"].ravel(), data["slope"].ravel(),
    ])
    matched = ci.covariance.covariance_halo_response(
        inputs=inputs, growth_coefficient=47.0/21.0,
        dilation_coefficient=1.0/3.0, fractional=False)[1].reshape(shape)
    reference = ((47.0/21.0 - data["slope"]/3) * data["i11"]**2
                 * data["linear"] + data["i12"])
    scale = np.abs(data["response_onecov"])
    if np.any(scale == 0):
        raise ValueError("zero response needs a different fractional-error metric")
    difference = np.abs(matched - data["response_onecov"]) / scale
    numpy_error = np.abs(reference - data["response_onecov"]) / scale

    inputs[5] = data["two_halo_slope"].ravel()
    two_halo = ci.covariance.covariance_halo_response(
        inputs=inputs, growth_coefficient=47.0/21.0,
        dilation_coefficient=1.0/3.0, fractional=False)[1].reshape(shape)
    selected = (data["k_h_Mpc"] >= 0.001) & (data["k_h_Mpc"] <= 10.0)
    rows = []
    for row, redshift in enumerate(data["redshift"]):
        rows.append({
            "redshift": float(redshift),
            "matched_kernel_max_fractional_difference": float(difference[row].max()),
            "numpy_max_fractional_difference": float(numpy_error[row].max()),
            "two_halo_slope_max_fractional_change_k_0p001_to_10": float(
                np.max(np.abs(two_halo[row, selected] / matched[row, selected] - 1))),
        })
    record = {
        "scope": "shared-ingredient matter-response formula, not projected SSC",
        "export": metadata, "rows": rows,
        "interface_sha256": sha256(filename=ci.__file__),
        "script_sha256": sha256(filename=__file__),
        "passed": bool(max(difference.max(), numpy_error.max()) < 1.e-11),
        "remaining": ["independent halo inputs", "fractional response transfer",
                      "survey window and radial projection", "SSC convergence"],
    }
    np.savez_compressed(file=output / "responses.npz", matched=matched,
                        onecov=data["response_onecov"], two_halo_slope=two_halo,
                        redshift=data["redshift"], k_h_Mpc=data["k_h_Mpc"])
    (output / "comparison.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(rows, indent=2))
    if not record["passed"]:
        raise SystemExit("Shared-ingredient response mismatch; results preserved")


def main():
    """Run one explicit stage in its own numerical environment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("export", "compare"))
    parser.add_argument("source", type=Path,
                        help="OneCov INI for export, ingredient folder for compare")
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output directory")
    if int(os.environ.get("OMP_NUM_THREADS", "0")) < 1:
        parser.error("set OMP_NUM_THREADS explicitly")
    args.output.mkdir(parents=True)
    if args.mode == "export":
        export_response(source=args.source.resolve(), output=args.output.resolve(),
                        onecov=args.onecov.resolve())
    else:
        compare_response(source=args.source.resolve(), output=args.output.resolve())


if __name__ == "__main__":
    main()
