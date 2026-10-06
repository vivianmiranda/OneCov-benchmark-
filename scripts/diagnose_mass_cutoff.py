"""Measure CoCoA's low-mass completion with its halo model unchanged.

Use the isolated wider-domain build for the cutoff scan, or the installed
interface for a baseline. A common sigma table and retained mass panels
keep the cutoff comparison separate from mass-table regridding. Outputs
are ingredients at three redshifts, not a full survey covariance.
"""

import argparse
import contextlib
import json
import math
import os
from pathlib import Path
import signal
import sys
import time

import numpy as np

from common import revision, sha256


def mass_edges(exponent, tail_panels=None):
    """Append low-mass panels without moving the original eight panels."""
    if tail_panels is not None and exponent < 4:
        lower = [exponent]
        if tail_panels == "intervals":
            lower = range(exponent, 4, 4)
        return np.concatenate(([math.log(10.0**value) for value in lower],
                               mass_edges(4)))
    original = np.linspace(np.log(1e6), np.log(1e17), 9)
    lower = np.log(10.0)*np.arange(exponent, 6)
    if len(lower):
        # Match the scalar C log used by the exact mass-domain guard.
        # Multiplying a negative exponent by log(10) can round below it.
        lower[0] = math.log(10.0**exponent)
    return np.concatenate((lower, original))


def completion(ci, scale, wave, edges, nquad, tail_nodes=None):
    """Read the missing weight from the same GSL rule and core halo fits.

    At zero wavenumber the normalized profile is one. Each mass node's
    contribution is dlnM*f(nu)*nu*dlnnu/dlnM*b(nu). Summing those weights
    measures how much of the required unit response is explicitly resolved.
    This diagnostic does not replace the production moment calculation.
    """
    nodes, weights = ci.covariance.covariance_integration_rule(nquad=nquad)
    masses = []
    measures = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        panel_nodes, panel_weights = nodes, weights
        if tail_nodes is not None and lower < math.log(1e4):
            # The private C build uses GSL's 32-point table when requested.
            # The public binding deliberately retains its >=64 contract.
            # Gauss-Legendre roots provide the same rule for this scalar
            # diagnostic of F and B; halo moments still come from C.
            panel_nodes, panel_weights = np.polynomial.legendre.leggauss(
                tail_nodes)
        half = (upper-lower)/2
        masses.extend(np.exp((upper+lower)/2+half*panel_nodes))
        measures.extend(half*panel_weights)
    masses = np.asarray(masses)
    measures = np.asarray(measures)
    missing = []
    resolved_mass = []
    profiles = []
    minimum = float(np.exp(edges[0]))
    for a in scale:
        contributions = []
        mass_contributions = []
        for mass, measure in zip(masses, measures):
            peak = 1.686/np.sqrt(ci.sigma2(float(mass), float(a), 1))
            slope = ci.dlognudlogm(float(mass), float(a))
            multiplicity = ci.fnu(float(peak), float(a))
            bias = ci.hb1nu(float(peak), float(a))
            contributions.append(measure*multiplicity*peak*slope*bias)
            mass_contributions.append(measure*multiplicity*peak*slope)
        missing.append(1.0-sum(contributions))
        resolved_mass.append(sum(mass_contributions))
        concentration = ci.conc(minimum, float(a))
        profiles.append([ci.u_nfw_c(concentration, float(k), minimum,
                                    float(a)) for k in wave])
    return (np.asarray(missing), np.asarray(profiles),
            np.asarray(resolved_mass))


def run(args):
    root = Path(os.environ["ROOTDIR"])
    core = root/"external_modules/code/cosmolike_core"
    sys.path[:0] = [str(core), str(root/"projects/lsst_y1/covariance")]
    build_record = None
    if args.interface:
        sys.path.insert(0, str(args.interface.resolve()))
        build_record = json.loads((args.interface/"build.json").read_text())
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.halo import halo_power_response
    from cosmolike_notebook_utils.covariance.geometry import angular_rule

    if args.interface and sha256(ci.__file__) != build_record["interface_sha256"]:
        raise ValueError("imported interface differs from the selected build")
    if args.tail_nodes is not None:
        if build_record is None or not build_record.get("split_tail"):
            raise ValueError("the split rule needs the isolated split-tail build")
        os.environ["COCOA_DIAGNOSTIC_TAIL_NQUAD"] = str(args.tail_nodes)
    # A private build can reproduce the old 1e6 table or extend it.
    # Do not request a mass integral below that build's supported domain.
    exponents = [6]
    if build_record is not None:
        for exponent in (4, 2):
            if 10.0**exponent >= build_record["table_mass_min"]:
                exponents.append(exponent)
    if args.exponents is not None:
        exponents = args.exponents
        supported_min = 1e4
        if build_record is not None:
            supported_min = build_record["table_mass_min"]
        if 10.0**min(exponents) < supported_min:
            raise ValueError("mass cutoff is below the selected table domain")
    settings = configuration(gaussian={"nonlimber": False, "ia": "none"},
                             integration_accuracy=2)
    tables = initialize(ci, settings)
    np.savez_compressed(args.output/"power_tables.npz", **tables)
    ci.init_accuracy_boost(accuracy_boost=args.boost, integration_accuracy=2)

    # Include both covariance diagonals and all unequal pairs. k=0 is
    # used only for the separate normalization check, not tree kernels.
    redshift = np.array([0.1, 0.5, 1.0])
    scale = 1/(1+redshift)
    length = 2997.92458
    wave_h = np.geomspace(1e-3, 100.0, 41)
    wave = wave_h*length
    grids = np.tile(wave, (len(scale), 1))
    first, second = np.triu_indices(len(wave))
    pairs = np.array([wave[first], wave[second]])
    _, angle_weight, corner = angular_rule(
        nquad=256, npanel=8, interface=ci.covariance)

    # Angular averages depend on power, not the mass cutoff. Compute
    # them once and supply exactly the same arrays to every cutoff.
    linear = []
    tree = []
    for a in scale:
        power = ci.covariance.covariance_power(a=float(a), k=wave, linear=True)
        pk = np.array([power[first], power[second]])
        magnitude = np.sqrt((pairs[0, :, None]-pairs[1, :, None])**2
                           +2*pairs[0, :, None]*pairs[1, :, None]*corner)
        internal = ci.covariance.covariance_power(
            a=float(a), k=np.ascontiguousarray(magnitude), linear=True)
        tree.append(ci.covariance.covariance_tree_averages(
            k=pairs, pk=pk, corner=corner, weight=angle_weight, ps=internal))
        linear.append(pk)

    saved = dict(redshift=redshift, k=wave_h, first=first, second=second)
    records = []
    for exponent in exponents:
        edges = mass_edges(exponent, args.tail_panels)
        for nquad in args.nodes:
            key = f"m{exponent}_n{nquad}"
            arguments = dict(a=scale, k=grids, lnm_edges=edges, nquad=nquad)
            start = time.perf_counter()
            single, moments = ci.covariance.covariance_halo_moments(**arguments)
            first_seconds = time.perf_counter()-start
            timings = []
            for _ in range(args.repeats):
                start = time.perf_counter()
                repeated = ci.covariance.covariance_halo_moments(**arguments)
                timings.append(time.perf_counter()-start)
                if not np.array_equal(single, repeated[0]):
                    raise ValueError("I11 is not repeatable")
                if not np.array_equal(moments, repeated[1]):
                    raise ValueError("pair moments are not repeatable")

            missing, profiles, resolved_mass = completion(
                ci, scale, wave, edges, nquad, args.tail_nodes)
            zero, _ = ci.covariance.covariance_halo_moments(
                a=scale, k=np.zeros((len(scale), 1)), lnm_edges=edges,
                nquad=nquad, pair_moments=False)
            if np.max(np.abs(zero-1)) > 1e-12:
                raise ValueError("the corrected large-scale limit failed")
            terms = []
            for row in range(len(scale)):
                terms.append(ci.covariance.covariance_halo_trispectrum(
                    pk=linear[row], tree=tree[row],
                    i11=np.array([single[row, first], single[row, second]]),
                    moments=np.ascontiguousarray(moments[:, row, :])))

            # The native response helper chooses the same GSL size via
            # integration_accuracy; its derivative step is held fixed.
            level = {96: 0, 128: 1, 256: 2, 512: 3}[nquad]
            response = halo_power_response(
                interface=ci.covariance, a=scale, k=grids, lnm_edges=edges,
                accuracy_boost=1, mnu=0, integration_accuracy=level)
            saved[key+"_i11"] = single
            saved[key+"_moments"] = moments
            saved[key+"_terms"] = np.asarray(terms)
            saved[key+"_response"] = response
            saved[key+"_missing"] = missing
            saved[key+"_correction"] = missing[:, None]*profiles
            records.append(dict(
                key=key, mass_min=10.0**exponent, mass_panels=len(edges)-1,
                nquad=nquad, missing_weight=missing.tolist(),
                resolved_mass=resolved_mass.tolist(),
                zero_mode_error=float(np.max(np.abs(zero-1))),
                first_call_seconds=first_seconds, seconds=timings,
            ))
            print(f"{key}: A={missing}; mean={np.mean(timings):.6f}s",
                  flush=True)

    mass = np.geomspace(10.0**min(exponents), 1e17, 301)
    mass = np.unique(np.concatenate((mass, 10.0**np.asarray(exponents))))
    saved["mass"] = mass
    saved["sigma"] = np.array([
        [np.sqrt(ci.sigma2(float(m), float(a), 1)) for m in mass]
        for a in scale])
    # Archive the same edge continuation for a separate native OneCov
    # filter check of the FFTLog integration boundary. The power model is
    # unchanged; extending an extrapolation does not calibrate that model.
    tail_decades = args.tail_log10_max+7
    tail_wave = np.geomspace(1e-7, 10.0**args.tail_log10_max,
                            tail_decades*4096+1)
    saved["tail_k"] = tail_wave
    saved["tail_power"] = np.array([
        ci.covariance.covariance_power(
            a=float(a), k=tail_wave*length, linear=True)*length**3
        for a in scale])
    saved["rho"] = 7.4775e21*settings["cosmology"]["omegam"]/length**3
    np.savez_compressed(args.output/"ingredients.npz", **saved)
    record = dict(
        core=revision(core), interface_sha256=sha256(ci.__file__),
        diagnostic_build=build_record, boost=args.boost, records=records,
        threads=int(os.environ["OMP_NUM_THREADS"]),
        power_table_sha256=sha256(args.output/"power_tables.npz"),
        ingredients_sha256=sha256(args.output/"ingredients.npz"),
        script_sha256=sha256(__file__),
        supplied_k_hmpc=[float(10**tables["log10k_2D"][0]),
                        float(10**tables["log10k_2D"][-1])],
        cosmology=settings["cosmology"],
        tail_nodes=args.tail_nodes, tail_panels=args.tail_panels,
    )
    (args.output/"report.json").write_text(json.dumps(record, indent=2)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interface", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--boost", type=float, default=4)
    parser.add_argument("--nodes", type=int, nargs="+", default=[96, 128, 256])
    parser.add_argument("--repeats", type=int, default=11)
    parser.add_argument("--exponents", type=int, nargs="+",
                        choices=(-20, -16, -12, -8, -4, -3, 0, 2, 4, 6))
    parser.add_argument("--tail-nodes", type=int, choices=(32, 64, 96, 128, 256))
    parser.add_argument("--tail-panels", choices=("single", "intervals"))
    parser.add_argument("--tail-log10-max", type=int, choices=(7, 9, 15), default=7)
    args = parser.parse_args()
    if args.output.exists() or args.boost not in (1, 2, 4, 8):
        parser.error("use a new output directory and boost 1,2,4 or 8")
    if not set(args.nodes) <= {96, 128, 256, 512} or args.repeats < 2:
        parser.error("use supported quadratures and at least two repeats")
    if (args.tail_nodes is None) != (args.tail_panels is None):
        parser.error("specify both tail-nodes and tail-panels, or neither")
    signal.alarm(600)
    args.output.mkdir(parents=True)
    with (args.output/"run.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            run(args)
    print(f"Saved cutoff diagnostic: {args.output}")


if __name__ == "__main__":
    main()
