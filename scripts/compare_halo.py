"""Compare native halo ingredients at z=0.1, 0.5 and 1.

The two codes retain their own halo fits and concentration relations.
Differences are physical/numerical diagnostics, not an equality test.
Mass is Msun/h, k is h/Mpc, number density is (h/Mpc)^3, and P is (Mpc/h)^3.
Export in the OneCov environment, then compare in Cocoa's environment.
"""

import argparse
import configparser
import contextlib
import json
import os
import signal
import sys
from pathlib import Path

import numpy as np

from common import power_table_record, revision, sha256


def export(args):
    """Read native halo statistics and moments without changing OneCov."""
    sys.path.insert(0, str(args.onecov.resolve()))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_polyspectra import PolySpectra

    config = configparser.ConfigParser()
    config.read(args.source)
    config["output settings"]["directory"] = str(args.output.resolve())
    config["halomodel evaluation"]["M_bins"] = str(args.mass_nodes)
    config["misc"]["num_cores"] = os.environ["OMP_NUM_THREADS"]
    filename = args.output/"onecov.ini"
    with filename.open("w") as stream:
        config.write(stream)
    rows = []
    redshifts = [0.1, 0.5, 1.0]
    with (args.output/"native.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
                Input().read_input(config_name=str(filename)))
            tables = FileInput(bias).read_input(config_name=str(filename))
            model = PolySpectra(redshifts[0], terms, obs, cosmo, bias, hod,
                                survey, prec, tables)
            for redshift in redshifts:
                model.update_mass_func(redshift, bias, hod, prec)
                mass = model.mass_func.m.copy()
                wave = model.mass_func.k.copy()
                peak = np.sqrt(model.mass_func.nu)
                row = dict(
                    mass=mass, k=wave, sigma=model.mass_func.sigma.copy(),
                    dndlnm=model.mass_func.dndm*mass,
                    bias=model.bias(bias, prec["hm"]).copy(),
                    raw_bias=model._HaloModel__bias_tinker10_fittfunc(peak),
                    concentration=model._HaloModel__concentration("duffy08"),
                    profile=model.uk(bias).copy(),
                    linear=model.mass_func.power.copy(),
                    rho=float(model.rho_bg), norm_bias=float(model.norm_bias))
                row["i11"] = model.halo_model_integral_I_alpha_x(
                    bias, hod, prec["hm"], 1, "m")[:, 0]
                row["i02"] = model._PolySpectra__P_xy_1h(
                    bias, hod, prec["hm"], "m", "m")[:, 0]
                row["i12"] = np.diagonal(model.halo_model_integral_I_alpha_xy(
                    bias, hod, prec["hm"], 1, "m", "m")[:, :, 0, 0]).copy()
                row["response"] = model.powspec_responses(
                    bias, hod, prec["hm"])[2][:, 0]
                rows.append(row)
    arrays = {name: np.array([row[name] for row in rows]) for name in rows[0]}
    np.savez_compressed(args.output/"inputs.npz", redshift=redshifts, **arrays)
    record = dict(onecov=revision(args.onecov), source_sha256=sha256(args.source),
                  config_sha256=sha256(filename), script_sha256=sha256(__file__),
                  inputs_sha256=sha256(args.output/"inputs.npz"),
                  mass_nodes=args.mass_nodes, redshifts=redshifts,
                  units="Msun/h, h/Mpc, (Mpc/h)^3",
                  native_choices="Tinker10, normalized halo bias, Duffy08 M200m; matter moment lower mass 1e2")
    (args.output/"export.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"Exported native halo ingredients: {args.output}")


def compare(args):
    """Evaluate Cocoa's public halo readers at the exported masses and k."""
    cocoa = args.cocoa.resolve()
    core = cocoa/"external_modules/code/cosmolike_core"
    sys.path[:0] = [str(core), str(cocoa/"projects/lsst_y1/covariance")]
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.halo import halo_power_response

    metadata = json.loads((args.source/"export.json").read_text())
    if sha256(args.source/"inputs.npz") != metadata["inputs_sha256"]:
        raise ValueError("halo inputs changed after export")
    data = np.load(args.source/"inputs.npz")
    settings = configuration(gaussian={"nonlimber": False, "ia": "none"},
                             integration_accuracy=args.integration_accuracy)
    tables = initialize(ci, settings)
    # Unit constants come from the initialized core's structs.c.
    length = 2997.92458
    rho = 7.4775e21*settings["cosmology"]["omegam"]/length**3
    values = {name: [] for name in ("sigma", "dndlnm", "bias", "concentration",
                                    "linear", "profile_matched_c", "bias_matched_nu")}
    scale = 1/(1+data["redshift"])
    for row, a in enumerate(scale):
        mass, wave = data["mass"][row], data["k"][row]
        sigma = np.sqrt([ci.sigma2(float(m), float(a), 1) for m in mass])
        peak = 1.686/sigma
        slope = np.array([ci.dlognudlogm(float(m), float(a)) for m in mass])
        multiplicity = np.array([ci.fnu(float(nu), float(a)) for nu in peak])
        values["sigma"].append(sigma)
        values["dndlnm"].append(rho/mass*multiplicity*peak*slope)
        values["bias"].append([ci.hb1nu(float(nu), float(a)) for nu in peak])
        values["concentration"].append([ci.conc(float(m), float(a)) for m in mass])
        values["linear"].append(ci.covariance.covariance_power(
            a=float(a), k=np.ascontiguousarray(wave*length), linear=True)*length**3)

        # Equal peak height and concentration isolate the bias formula and
        # NFW transform. They are explicitly separate from native choices.
        values["bias_matched_nu"].append([
            ci.hb1nu(float(nu), float(a)) for nu in 1.686/data["sigma"][row]])
        values["profile_matched_c"].append([
            [ci.u_nfw_c(float(c), float(k*length), float(m), float(a))
             for m, c in zip(mass, data["concentration"][row])]
            for k in wave])
    values = {name: np.asarray(value) for name, value in values.items()}

    # Each batch row contains just one k. This requests diagonal I02/I12
    # without spending time on pair moments not used by this comparison.
    repeated_a = np.repeat(scale, data["k"].shape[1])
    single, moments = ci.covariance.covariance_halo_moments(
        a=repeated_a, k=np.ascontiguousarray(data["k"].ravel()[:, None]*length),
        lnm_edges=settings["lnm_edges"], nquad=settings["halo_mass_nquad"])
    shape = data["k"].shape
    values["i11"] = single.reshape(shape)
    values["i02"] = moments[0].reshape(shape)*length**3
    values["i12"] = moments[1].reshape(shape)*length**3
    values["response"] = halo_power_response(
        interface=ci.covariance, a=scale, k=np.ascontiguousarray(data["k"]*length),
        lnm_edges=settings["lnm_edges"], accuracy_boost=1, mnu=0,
        integration_accuracy=args.integration_accuracy)*length**3
    np.savez_compressed(args.output/"cocoa.npz", **values)

    rows = []
    for row, redshift in enumerate(data["redshift"]):
        mass = (data["mass"][row] >= 1.e10) & (data["mass"][row] <= 1.e15)
        wave = (data["k"][row] >= 0.001) & (data["k"][row] <= 10)
        report = {"redshift": float(redshift)}
        for name in ("sigma", "dndlnm", "bias", "concentration", "linear",
                     "i11", "i02", "i12", "response"):
            chosen = mass if name in ("sigma", "dndlnm", "bias", "concentration") else wave
            difference = values[name][row, chosen]/data[name][row, chosen]-1
            report[name+"_fractional_range"] = [float(difference.min()), float(difference.max())]
        report["raw_bias_matched_peak_max_fractional"] = float(np.max(np.abs(
            values["bias_matched_nu"][row]/data["raw_bias"][row]-1)))
        # An absolute profile error avoids magnifying oscillatory near-zeros.
        report["nfw_matched_concentration_max_absolute"] = float(np.max(np.abs(
            values["profile_matched_c"][row]-data["profile"][row])))
        rows.append(report)
    record = dict(scope="native halo prescriptions; equality is not expected",
                  installed_power=power_table_record(tables, settings),
                  export=metadata, rows=rows, integration_accuracy=args.integration_accuracy,
                  core=revision(core), interface_sha256=sha256(ci.__file__),
                  script_sha256=sha256(__file__),
                  statistic_mass_range=[1.e10, 1.e15], statistic_k_range=[0.001, 10],
                  cocoa_density=rho, onecov_density=data["rho"].tolist(),
                  mass_min=float(np.exp(settings["lnm_edges"][0])),
                  mass_max=float(np.exp(settings["lnm_edges"][-1])),
                  native_choices="Tinker10 bias-consistent multiplicity, Bhattacharya13, M200m; production mass panels and I11 completion; fractional response transferred to nonlinear P")
    (args.output/"comparison.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(rows, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("export", "compare"))
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--cocoa", type=Path,
                        default=Path(__file__).resolve().parents[2]/"cocoa/Cocoa")
    parser.add_argument("--mass-nodes", type=int, default=200)
    parser.add_argument("--integration-accuracy", type=int, choices=range(5), default=0)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output directory")
    args.output.mkdir(parents=True)
    signal.alarm(args.timeout)
    if args.mode == "export":
        export(args)
    else:
        compare(args)


if __name__ == "__main__":
    main()
