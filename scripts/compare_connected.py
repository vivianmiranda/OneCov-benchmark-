"""Compare a shear cNG projection using OneCov's actual matter tables.

Eight multipoles between 30 and 3000 test the radial projection, including
off-diagonal entries. The native five-band G+cNG matrix is saved separately
for positivity and refinement checks. Distances are Mpc/h, T is (Mpc/h)^9.
No covariance source is changed and no native-model equality is assumed.
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
    """Compute native G+cNG and capture the cNG radial integrand."""
    sys.path.insert(0, str(args.onecov.resolve()))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_ell_space import CovELLSpace

    config = configparser.ConfigParser()
    config.read(args.source)
    config["output settings"]["directory"] = str(args.output.resolve())
    config["covariance terms"].update(nongauss="True", ssc="False")
    config["misc"]["num_cores"] = os.environ["OMP_NUM_THREADS"]
    config["covELLspace settings"].update(
        ell_min="30", ell_max="3000", ell_bins="100",
        integration_steps=str(args.radial_nodes), tri_delta_z=str(args.delta_z))
    config["halomodel evaluation"]["M_bins"] = str(args.mass_nodes)
    config["trispec evaluation"].update(
        log10k_bins=str(args.k_nodes), log10k_min="-4", log10k_max="2",
        matter_klim=str(args.corner_guard), matter_mulim=str(args.corner_guard))
    filename = args.output/"onecov.ini"
    with filename.open("w") as stream:
        config.write(stream)

    started = time.perf_counter()
    with (args.output/"native.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
                Input().read_input(config_name=str(filename)))
            tables = FileInput(bias).read_input(config_name=str(filename))
            model = CovELLSpace(terms, obs, output, cosmo, bias, ia, hod,
                               survey, prec, tables)
            if model.gg or model.gm or model.n_tomo_lens != 1:
                raise ValueError("select one source population without lenses")
            (gauss, connected, _), captured = capture_return(
                model.covELL_non_gaussian,
                ("trispec_integrand_mmmm", "nongaussELLmmmm", "zet_list"),
                lambda: model.calc_covELL(obs, output, bias, hod, survey,
                                         prec, tables))

    chi = model.los_integration_chi
    chosen = np.linspace(0, len(model.ellrange)-1, 8, dtype=int)
    area = float(survey["survey_area_lens"][0])*(np.pi/180)**2
    # The public full calculation adds area normalization during binning.
    # Its captured unbinned cNG has not yet been divided by survey area.
    raw = captured["nongaussELLmmmm"][:, :, 0, 0, 0, 0, 0, 0]/area
    native = raw[np.ix_(chosen, chosen)]
    trispectrum = captured["trispec_integrand_mmmm"][:, :, :, 0, 0]
    trispectrum = trispectrum[:, chosen][:, :, chosen]
    band_g = sum(gauss[15:18])[:, :, 0, 0, 0, 0, 0, 0]
    band_cng = connected[5][:, :, 0, 0, 0, 0, 0, 0]
    target = args.output/"inputs.npz"
    np.savez_compressed(
        target, chi=chi, dchi=simpson(np.eye(len(chi)), x=chi, axis=1),
        window=model.spline_lensweight[0](chi), ell=model.ellrange[chosen],
        trispectrum=trispectrum, onecov=native, area_sr=area,
        band_ell=model.ellrange_lensing, band_gaussian=band_g,
        band_connected=band_cng, band_total=band_g+band_cng)
    total = band_g+band_cng
    normalized = total/np.sqrt(np.outer(total.diagonal(), total.diagonal()))
    eigenvalues = np.linalg.eigvalsh(normalized)
    if not np.isfinite(total).all() or eigenvalues[0] <= 0:
        raise ValueError("native five-band G+cNG is not finite positive definite")
    record = dict(onecov=revision(args.onecov), inputs_sha256=sha256(target),
                  script_sha256=sha256(__file__), config_sha256=sha256(filename),
                  observer_sha256=sha256(Path(__file__).with_name("capture_native.py")),
                  radial_nodes=len(chi), k_nodes=args.k_nodes,
                  mass_nodes=args.mass_nodes, delta_z=args.delta_z,
                  corner_guard=args.corner_guard,
                  redshifts=captured["zet_list"],
                  total_correlation_eigenvalues=eigenvalues.tolist(),
                  instrumented_seconds=time.perf_counter()-started,
                  timing_note="diagnostic capture, not an uninstrumented benchmark")
    (args.output/"export.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))


def compare(args):
    """Project shared T(k,q,chi) with Cocoa's production C assembly."""
    import cosmolike_lsst_y1_interface as ci

    metadata = json.loads((args.source/"export.json").read_text())
    if sha256(args.source/"inputs.npz") != metadata["inputs_sha256"]:
        raise ValueError("native export changed")
    data = np.load(args.source/"inputs.npz")
    size = len(data["ell"])
    projected = np.zeros((4*size, 4*size, len(data["chi"])))
    # These are Fourier samples, so no spin transform remains to apply.
    # Use one of the four generic operator slots and leave the others zero.
    projected[:size, :size] = data["trispectrum"].transpose(1, 2, 0)
    result = ci.covariance.covariance_project_connected(
        probes=np.array([0], dtype=np.int32), projected=projected,
        pair_window=np.ascontiguousarray(data["window"][None, :]**2),
        measure=data["dchi"]/(data["area_sr"]*data["chi"]**6))
    reference = data["onecov"]
    normalized = (result-reference)/np.sqrt(np.outer(
        reference.diagonal(), reference.diagonal()))
    if not np.isfinite(result).all():
        raise ValueError("nonfinite Cocoa projection")
    np.savez_compressed(args.output/"comparison.npz", cocoa=result,
                        onecov=reference, ell=data["ell"])
    record = dict(export=metadata, scope="shared-input radial cNG projection",
                  max_variance_scaled_residual=float(np.max(np.abs(normalized))),
                  max_fractional_residual=float(np.max(np.abs(result/reference-1))),
                  symmetry_max_absolute=float(np.max(np.abs(result-result.T))),
                  interface_sha256=sha256(ci.__file__), script_sha256=sha256(__file__))
    (args.output/"comparison.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("export", "compare"))
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--k-nodes", type=int, default=9)
    parser.add_argument("--mass-nodes", type=int, default=400)
    parser.add_argument("--radial-nodes", type=int, default=300)
    parser.add_argument("--delta-z", type=float, default=0.5)
    parser.add_argument("--corner-guard", type=float, default=0.001)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output directory")
    args.output.mkdir(parents=True)
    signal.alarm(args.timeout)
    (export if args.mode == "export" else compare)(args)


if __name__ == "__main__":
    main()
