"""Compare separated matter trispectra, with native and shared halo inputs.

Export calls unchanged OneCov methods at three redshifts. Compare calls
CoCoA's production kernels. Units throughout are h/Mpc for k and
(Mpc/h)^9 for T. Upper-triangular pairs include unequal wavenumbers.
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
from scipy.interpolate import UnivariateSpline

from capture_native import capture_return
from common import power_table_record, revision, sha256


def export(args):
    """Read native terms and their actual inputs without replacing methods."""
    sys.path.insert(0, str(args.onecov.resolve()))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_polyspectra import PolySpectra

    config = configparser.ConfigParser()
    config.read(args.source)
    config["output settings"]["directory"] = str(args.output.resolve())
    config["covariance terms"].update(nongauss="True", ssc="False")
    config["halomodel evaluation"]["M_bins"] = str(args.mass_nodes)
    config["misc"]["num_cores"] = os.environ["OMP_NUM_THREADS"]
    config["trispec evaluation"].update(
        log10k_bins=str(args.k_nodes), log10k_min="-3", log10k_max="1",
        matter_klim=str(args.corner_guard),
        matter_mulim=str(args.corner_guard))
    filename = args.output/"onecov.ini"
    with filename.open("w") as stream:
        config.write(stream)

    redshifts = [0.1, 0.5, 1.0]
    rows = []
    started = time.perf_counter()
    with (args.output/"native.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
                Input().read_input(config_name=str(filename)))
            tables = FileInput(bias).read_input(config_name=str(filename))
            model = PolySpectra(redshifts[0], terms, obs, cosmo, bias, hod,
                                survey, prec, tables)
            for redshift in redshifts:
                model.update_mass_func(redshift, bias, hod, prec)
                _, values = capture_return(
                    model._PolySpectra__trispectra_234h,
                    ("integral_m", "integral_mm", "integral_mmm",
                     "Pspline_eval", "trispec_2h", "trispec_3h", "trispec_4h"),
                    lambda: model.trispectra(output, bias, hod,
                                            prec["hm"], tables["tri"]))
                values["onecov"] = np.array([
                    model.trispec1h_mmmm[:, :, 0, 0],
                    values.pop("trispec_2h"), values.pop("trispec_3h"),
                    values.pop("trispec_4h")])

                # OneCov stores these angle integrals divided by powers
                # of growth. Restore their values at this redshift. Its
                # first two integrals contain both +/- angles: divide by
                # two to obtain Cocoa's averages <P> and <B>.
                growth = model.mass_func.growth_factor
                values["tree"] = np.array([
                    model.int_2h*growth**2/2, model.int_3h*growth**4/2,
                    model.int_4h*growth**6])
                values["power_k"] = model.mass_func.k.copy()
                values["power"] = model.mass_func.power.copy()
                rows.append(values)
    arrays = {name: np.array([row[name] for row in rows]) for name in rows[0]}
    if not all(np.isfinite(value).all() for value in arrays.values()):
        raise ValueError("nonfinite native trispectrum ingredient")
    target = args.output/"inputs.npz"
    np.savez_compressed(target, redshift=redshifts, k=model.krange_tri, **arrays)
    record = dict(onecov=revision(args.onecov), inputs_sha256=sha256(target),
                  script_sha256=sha256(__file__), config_sha256=sha256(filename),
                  observer_sha256=sha256(Path(__file__).with_name("capture_native.py")),
                  mass_nodes=args.mass_nodes, k_nodes=len(model.krange_tri),
                  corner_guard=args.corner_guard, redshifts=redshifts,
                  instrumented_seconds=time.perf_counter()-started,
                  timing_note="diagnostic capture, not an uninstrumented benchmark")
    (args.output/"export.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))


def combine(terms):
    """Combine Cocoa's two distinct 2h partitions into the four halo orders."""
    return np.array([terms[0], terms[1]+terms[2], terms[3], terms[4]])


def compare(args):
    """Test shared-input assembly, angular averages, and native halo terms."""
    cocoa = args.cocoa.resolve()
    core = cocoa/"external_modules/code/cosmolike_core"
    sys.path[:0] = [str(core), str(cocoa/"projects/lsst_y1/covariance")]
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.geometry import angular_rule
    from cosmolike_notebook_utils.covariance.halo import halo_trispectrum

    metadata = json.loads((args.source/"export.json").read_text())
    if sha256(args.source/"inputs.npz") != metadata["inputs_sha256"]:
        raise ValueError("native export changed")
    data = np.load(args.source/"inputs.npz")
    settings = configuration(gaussian={"nonlimber": False, "ia": "none"},
                             integration_accuracy=args.integration_accuracy)
    tables = initialize(ci, settings)
    first, second = np.triu_indices(len(data["k"]))
    pairs = np.array([data["k"][first], data["k"][second]])
    _, weight, corner = angular_rule(nquad=settings["tree_nquad"],
                                     npanel=settings["tree_npanel"],
                                     interface=ci.covariance)
    saved = {name: [] for name in ("shared", "repeated_partition", "native",
                                   "tree", "onecov")}
    reports = []
    for row, redshift in enumerate(data["redshift"]):
        pk = data["Pspline_eval"][row][[first, second]]
        i11 = data["integral_m"][row][[first, second]]
        i13 = data["integral_mmm"][row]
        moments = np.array([
            np.zeros(len(first)), data["integral_mm"][row, first, second],
            i13[first, second], i13[second, first],
            data["onecov"][row, 0, first, second]])
        tree = np.ascontiguousarray(data["tree"][row][:, first, second])
        reference = data["onecov"][row][:, first, second]
        shared = combine(ci.covariance.covariance_halo_trispectrum(
            pk=pk, i11=i11, moments=moments, tree=tree))

        # Diagnose the two 1+3 halo partitions independently. This second
        # call repeats OneCov's use of I13(K,Q,Q) in both terms. It changes
        # only supplied diagnostic inputs, never either code's source.
        moments[3] = moments[2]
        repeated = combine(ci.covariance.covariance_halo_trispectrum(
            pk=pk, i11=i11, moments=moments, tree=tree))

        # Use exactly OneCov's quadratic log-power interpolation to supply
        # internal P(|K+Q|). Cocoa then performs its own angle integration.
        spline = UnivariateSpline(np.log(data["power_k"][row]),
                                   np.log(data["power"][row]), k=2, s=0)
        # Equal k pairs approach zero internal k. OneCov cuts out part of
        # that corner; feeding its extrapolated spline to Cocoa's finer
        # angular rule is not a shared-domain test. Compare unequal pairs
        # here, requiring every internal k to remain inside the input table.
        unequal = first != second
        selected = np.ascontiguousarray(pairs[:, unequal])
        magnitude = np.sqrt((selected[0]-selected[1])[:, None]**2
                           +2*(selected[0]*selected[1])[:, None]*corner)
        if (magnitude.min() < data["power_k"][row, 0]
                or magnitude.max() > data["power_k"][row, -1]):
            raise ValueError("angular comparison would extrapolate power")
        angular = ci.covariance.covariance_tree_averages(
            k=selected, pk=np.ascontiguousarray(pk[:, unequal]),
            corner=corner, weight=weight,
            ps=np.exp(spline(np.log(magnitude))))
        native = combine(halo_trispectrum(
            interface=ci.covariance, a=float(1/(1+redshift)),
            k=data["k"]*2997.92458, lnm_edges=settings["lnm_edges"],
            accuracy_boost=1, mnu=0,
            integration_accuracy=args.integration_accuracy)["terms"])*2997.92458**9
        for name, value in zip(saved, (shared, repeated, native, angular, reference)):
            if not np.isfinite(value).all():
                raise ValueError("nonfinite comparison term")
            saved[name].append(value)
        # Unequal grid pairs never sample below the supplied power range.
        # Equal pairs approach |K+Q|=0: OneCov's corner guard and quadratic
        # extrapolation then differ from Cocoa's production power reader.
        # Keep that diagnostic separate instead of calling it convergence.
        reports.append(dict(
            redshift=float(redshift),
            shared_max_fractional=np.max(np.abs(shared/reference-1), axis=1).tolist(),
            shared_diagonal_max_fractional=np.max(np.abs(
                shared[:, first == second]/reference[:, first == second]-1), axis=1).tolist(),
            repeated_partition_max_fractional=np.max(
                np.abs(repeated/reference-1), axis=1).tolist(),
            tree_offdiagonal_max_fractional=np.max(np.abs(
                angular/tree[:, unequal]-1), axis=1).tolist(),
            native_fractional_range=np.array([
                np.min(native/reference-1, axis=1),
                np.max(native/reference-1, axis=1)]).T.tolist()))
    np.savez_compressed(args.output/"comparison.npz", first=first, second=second,
                        k=data["k"], redshift=data["redshift"], **saved)
    record = dict(export=metadata, rows=reports, halo_order=["1h", "2h", "3h", "4h"],
                  installed_power=power_table_record(tables, settings),
                  integration_accuracy=args.integration_accuracy,
                  core=revision(core), interface_sha256=sha256(ci.__file__),
                  script_sha256=sha256(__file__))
    (args.output/"comparison.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(reports, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("export", "compare"))
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--cocoa", type=Path,
                        default=Path(__file__).resolve().parents[2]/"cocoa/Cocoa")
    parser.add_argument("--mass-nodes", type=int, default=400)
    parser.add_argument("--k-nodes", type=int, default=9)
    parser.add_argument("--corner-guard", type=float, default=0.001)
    parser.add_argument("--integration-accuracy", type=int, choices=range(5), default=0)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output directory")
    args.output.mkdir(parents=True)
    signal.alarm(args.timeout)
    (export if args.mode == "export" else compare)(args)


if __name__ == "__main__":
    main()
