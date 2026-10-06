"""Separate power, quadrature, density and halo-abundance differences.

Export Cocoa's actual linear reader and FFTLog variance in its environment.
Then call the unchanged hmf TopHat filter used by OneCov on these common
inputs. No independent variance implementation or numerical source patch
is introduced. Mass is Msun/h, radius Mpc/h and k h/Mpc throughout.
"""

import argparse
import contextlib
import json
import signal
import sys
from pathlib import Path

import numpy as np

from common import revision, sha256
from capture_native import capture_return


def cocoa(args, native):
    root = Path(__file__).resolve().parents[2]/"cocoa/Cocoa"
    core = root/"external_modules/code/cosmolike_core"
    sys.path[:0] = [str(core), str(root/"projects/lsst_y1/covariance")]
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.camb_cosmology import get_camb_cosmology

    settings = configuration(gaussian={"nonlimber": False, "ia": "none"},
                             integration_accuracy=2)
    # Observe the actual CAMB call used by the forecast. Reconstructing a
    # second call from a subset of its arguments can change the amplitude.
    tables, captured = capture_return(get_camb_cosmology, ("results", "pars"),
                                      lambda: initialize(ci, settings))
    selected = (native["mass"][0] >= 1e10) & (native["mass"][0] <= 1e15)
    mass = native["mass"][0, selected]
    redshift = native["redshift"]
    length = 2997.92458
    rho = 7.4775e21*settings["cosmology"]["omegam"]/length**3
    saved = dict(mass=mass, redshift=redshift, rho=rho)

    # Keep the supplied power fixed. This refinement changes only the
    # internal FFTLog output tables and their mass/time interpolation.
    for boost in (1, 2, 4):
        ci.init_accuracy_boost(accuracy_boost=float(boost), integration_accuracy=2)
        saved[f"sigma_boost{boost}"] = np.array([
            [np.sqrt(ci.sigma2(float(m), float(1/(1+z)), 1)) for m in mass]
            for z in redshift])
        if boost == 1:
            saved["slope"] = np.array([
                [ci.dlognudlogm(float(m), float(1/(1+z))) for m in mass]
                for z in redshift])
            nu_c = 1.686/saved["sigma_boost1"]
            nu_o = 1.686/native["sigma"][:, selected]
            for name, peak in (("fnu_native", nu_c), ("fnu_onecov_peak", nu_o)):
                saved[name] = np.array([
                    [ci.fnu(float(n), float(1/(1+z))) for n in peak[row]]
                    for row, z in enumerate(redshift)])

    # The native 200 nodes, a nested refinement of the same interval, and
    # a wide interval separate spacing errors from missing power tails.
    grids = dict(native=native["k"][0],
                 dense=np.geomspace(native["k"][0, 0], native["k"][0, -1],
                                    (native["k"].shape[1]-1)*64+1),
                 wide=np.geomspace(1e-7, 1e5, 65537))
    for name, wave in grids.items():
        saved[f"k_{name}"] = wave
        saved[f"power_{name}"] = np.array([
            ci.covariance.covariance_power(a=float(1/(1+z)),
                k=np.ascontiguousarray(wave*length), linear=True)*length**3
            for z in redshift])
    # CAMB's reported sigma8, its interpolated spectrum, and the core
    # reader are distinct numerical paths. Save the latter two at z=0
    # so OneCov's filter can test the amplitude conversion directly.
    wave = np.geomspace(1e-5, 100.0, 32769)
    saved["k_amplitude"] = wave
    # The covariance reader excludes a=1. Its nearest representable value
    # below one gives the z=0 limit without violating that public contract.
    saved["power_reader_z0"] = ci.covariance.covariance_power(
        a=float(np.nextafter(1.0, 0.0)),
        k=np.ascontiguousarray(wave*length), linear=True)*length**3
    camb_power = captured["results"].get_matter_power_interpolator(
        var1="delta_tot", var2="delta_tot", nonlinear=False,
        hubble_units=True, k_hunit=True, extrap_kmax=1000.0)
    saved["power_camb_z0"] = camb_power.P(0.0, wave)
    sigma8_fft = np.sqrt(ci.sigma2(float(4*np.pi*rho*8**3/3), 1.0, 1))
    np.savez_compressed(args.output/"inputs.npz", **saved)
    return dict(core=revision(core), interface_sha256=sha256(ci.__file__),
                cosmology=settings["cosmology"],
                camb_sigma8=float(captured["results"].get_sigma8_0()),
                core_sigma8_fft=float(sigma8_fft),
                camb_neff=float(captured["pars"].N_eff),
                camb_k_per_logint=int(captured["pars"].Transfer.k_per_logint),
                camb_kmax=float(captured["pars"].Transfer.kmax),
                power_table_shape=[len(tables["z_2D"]), len(tables["log10k_2D"])],
                power_table_sha256=sha256(args.output/"inputs.npz"))


def evaluate_filter(wave, power, radius):
    """Call OneCov's hmf filter in small batches to bound temporary memory."""
    from hmf.density_field.filters import TopHat
    filt = TopHat(wave, power)
    sigma = np.concatenate([filt.sigma(radius[start:start+32])
                            for start in range(0, len(radius), 32)])
    slope = np.concatenate([-0.5*filt.dlnss_dlnm(radius[start:start+32])
                            for start in range(0, len(radius), 32)])
    return sigma, slope


def onecov(args, native):
    root = Path(__file__).resolve().parents[2]/"OneCovariance"
    sys.path.insert(0, str(root))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_polyspectra import PolySpectra
    from hmf.density_field import filters

    cocoa_record = json.loads((args.cocoa/"report.json").read_text())
    if sha256(args.cocoa/"inputs.npz") != cocoa_record["power_table_sha256"]:
        raise ValueError("Cocoa inputs changed after export")
    shared = np.load(args.cocoa/"inputs.npz")
    selected = (native["mass"][0] >= 1e10) & (native["mass"][0] <= 1e15)
    mass = shared["mass"]
    radius_c = (3*mass/(4*np.pi*float(shared["rho"])))**(1/3)
    saved = {}
    rows = []
    for row, z in enumerate(native["redshift"]):
        radius_o = (3*mass/(4*np.pi*native["rho"][row]))**(1/3)
        sigma_o, slope_o = evaluate_filter(native["k"][row],
                                            native["linear"][row], radius_o)
        error = np.max(np.abs(sigma_o/native["sigma"][row, selected]-1))
        if error > 1e-12:
            raise ValueError("filter does not reproduce the native halo sigma")

        # Each replacement changes only one ingredient relative to the
        # preceding stage. Products of the ratios recover the native offset.
        stages = {"onecov_native": sigma_o}
        for name in ("native", "dense", "wide"):
            stages[f"cocoa_power_{name}_onecov_radius"] = evaluate_filter(
                shared[f"k_{name}"], shared[f"power_{name}"][row], radius_o)[0]
        for stride in (4, 2, 1):
            stages[f"shared_filter_stride{stride}"] = evaluate_filter(
                shared["k_wide"][::stride],
                shared["power_wide"][row, ::stride], radius_c)[0]
        for boost in (1, 2, 4):
            stages[f"cocoa_fft_boost{boost}"] = shared[f"sigma_boost{boost}"][row]
        for name, values in stages.items():
            saved.setdefault(name, []).append(values)

        # Factor dn/dlnM into density, f(nu), peak height and Jacobian.
        # Evaluating Cocoa f at OneCov's nu isolates its normalization.
        peak_o = 1.686/sigma_o
        peak_c = 1.686/shared["sigma_boost1"][row]
        f_o = native["dndlnm"][row, selected]*mass/(
            native["rho"][row]*peak_o*slope_o)
        factors = dict(density=np.full(len(mass), shared["rho"]/native["rho"][row]),
                       multiplicity_at_same_peak=shared["fnu_onecov_peak"][row]/f_o,
                       multiplicity_peak_shift=(shared["fnu_native"][row]
                                                /shared["fnu_onecov_peak"][row]),
                       peak_height=peak_c/peak_o,
                       mass_derivative=shared["slope"][row]/slope_o)
        abundance_ratio = np.prod(list(factors.values()), axis=0)
        for name, values in factors.items():
            saved.setdefault("abundance_"+name, []).append(values)
        saved.setdefault("abundance_ratio", []).append(abundance_ratio)
        rows.append(dict(z=float(z), native_filter_check=float(error),
                         factor_ranges={name: [float(v.min()), float(v.max())]
                                        for name, v in factors.items()},
                         abundance_ratio_range=[float(abundance_ratio.min()),
                                                float(abundance_ratio.max())]))

    # Refine OneCov's own power sampling, which also refines its numerical
    # sigma8 normalization. Halo prescriptions and cosmology stay fixed.
    terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
        Input().read_input(config_name=str(args.inputs/"onecov.ini")))
    tables = FileInput(bias).read_input(config_name=str(args.inputs/"onecov.ini"))
    model = PolySpectra(0.1, terms, obs, cosmo, bias, hod, survey, prec, tables)
    initial_step = model.mass_func.dlnk
    for refinement in (1, 2, 4, 8, 16, 32, 64):
        model.mass_func.update(dlnk=initial_step/refinement)
        sigmas = []
        slopes = []
        amplitudes = []
        sigma8 = []
        for row, z in enumerate(native["redshift"]):
            model.mass_func.update(z=float(z))
            radius = (3*mass/(4*np.pi*model.mass_func.mean_density0))**(1/3)
            sigma, slope = evaluate_filter(model.mass_func.k, model.mass_func.power, radius)
            sigmas.append(sigma)
            slopes.append(slope)
            amplitudes.append(float(model.mass_func._normalisation))
            sigma8.append(float(filters.TopHat(
                model.mass_func.k, model.mass_func.power).sigma(np.array([8.0]))[0]))
        saved[f"onecov_refine{refinement}"] = sigmas
        saved[f"onecov_slope_refine{refinement}"] = slopes
        saved[f"onecov_normalization_refine{refinement}"] = amplitudes
        saved[f"onecov_sigma8_refine{refinement}"] = sigma8

    saved = {name: np.asarray(values) for name, values in saved.items()}
    if np.max(np.abs(saved["onecov_refine1"]/saved["onecov_native"]-1)) > 1e-12:
        raise ValueError("refinement baseline differs from archived OneCov")
    np.savez_compressed(args.output/"diagnostic.npz", mass=mass, **saved)
    return dict(onecov=revision(root), rows=rows,
                growth_model=type(model.mass_func.growth).__name__,
                sigma8_input=float(cosmo["sigma8"]),
                growth_factors=[float(model.mass_func.growth.growth_factor(float(z)))
                                for z in native["redshift"]],
                filter_source_sha256=sha256(filters.__file__),
                cocoa_export=cocoa_record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("cocoa", "onecov"))
    parser.add_argument("--inputs", type=Path, default=Path("work/halo_onecov_800"))
    parser.add_argument("--cocoa", type=Path, default=Path("work/sigma_cocoa"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.output.exists() or args.timeout <= 0:
        parser.error("choose a new output directory and a positive timeout")
    signal.alarm(args.timeout)
    args.output.mkdir(parents=True)
    meta = json.loads((args.inputs/"export.json").read_text())
    if sha256(args.inputs/"inputs.npz") != meta["inputs_sha256"]:
        raise ValueError("native inputs changed")
    native = np.load(args.inputs/"inputs.npz")
    with (args.output/"run.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            record = cocoa(args, native) if args.mode == "cocoa" else onecov(args, native)
    filename = "inputs.npz" if args.mode == "cocoa" else "diagnostic.npz"
    record.update(native_export=meta, script_sha256=sha256(__file__),
                  output_sha256=sha256(args.output/filename))
    (args.output/"report.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"Saved {args.mode} sigma diagnostics to {args.output}")


if __name__ == "__main__":
    main()
