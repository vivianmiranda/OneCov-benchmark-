"""Compare every entry and variance mode of two complete saved covariances.

Gaussian, SSC, cNG and total differences remain separate. A component may
be singular or signed: its difference is measured relative to the positive
total reference, rather than attempting to invert the component itself.
"""

import argparse
import json
import os
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import cholesky, eigvalsh
from threadpoolctl import threadpool_limits

from common import sha256


def read_covariance(folder):
    """Require a complete archived run with its recorded input fingerprints."""
    record = json.loads((folder/"report.json").read_text())
    filename = folder/"covariance.npz"
    if sha256(filename) != record["covariance_sha256"]:
        raise ValueError(f"covariance changed after generation: {folder}")
    return dict(np.load(filename, allow_pickle=False)), record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--label", help="describe a recorded model diagnostic")
    args = parser.parse_args()
    if args.output.exists() or args.figure.exists():
        parser.error("choose new output and figure paths")
    candidate, first_record = read_covariance(args.candidate)
    reference, second_record = read_covariance(args.reference)
    if first_record["power_table_sha256"] != second_record["power_table_sha256"]:
        raise ValueError("power-spectrum inputs differ between the runs")
    for name in ("rows", "coordinate", "geometry", "coarse_ell",
                 "pair_area_sr2"):
        np.testing.assert_array_equal(candidate[name], reference[name])
    first_settings = json.loads(str(candidate["settings_json"]))
    second_settings = json.loads(str(reference["settings_json"]))
    # Only the mass panels and the recorded execution diagnostic may
    # change. In particular, angular/radial sampling and cosmology must not.
    controls = []
    for settings in (first_settings, second_settings):
        controls.append({key: value for key, value in settings.items()
                         if key not in ("lnm_edges", "execution")})
    if controls[0] != controls[1]:
        raise ValueError("settings other than the mass cutoff/domain differ")
    names = ("gaussian", "ssc", "cng", "total")
    diagonal = np.diag(reference["total"])
    if np.any(diagonal <= 0):
        raise ValueError("reference total has a nonpositive diagonal")
    normalization = np.sqrt(diagonal[:, None]*diagonal[None, :])
    for data in (candidate, reference):
        for name in names:
            values = data[name]
            if (values.shape != normalization.shape
                    or not np.all(np.isfinite(values))):
                raise ValueError(f"invalid full {name} matrix")
            if np.max(np.abs(values-values.T)/normalization) > 1e-12:
                raise ValueError(f"asymmetric full {name} matrix")
        summed = data["gaussian"]+data["ssc"]+data["cng"]
        if np.max(np.abs(summed-data["total"])/normalization) > 1e-12:
            raise ValueError("components do not sum to the archived total")

    report = dict(
        ndata=len(diagonal), masked_entries=0,
        candidate=first_record, reference=second_record, components={},
        candidate_settings=first_settings, reference_settings=second_settings,
        script_sha256=sha256(__file__),
    )
    normalized_reference = reference["total"]/normalization
    normalized_candidate = candidate["total"]/normalization
    with threadpool_limits(limits=1, user_api="blas"):
        # Cholesky checks positivity without altering eigenvalues. The
        # diagonal rescaling improves conditioning but preserves the sign.
        cholesky(normalized_reference, lower=True)
        cholesky(normalized_candidate, lower=True)
        report["both_positive_definite"] = True
        report["minimum_normalized_eigenvalues"] = [
            float(eigvalsh(matrix, subset_by_index=[0, 0])[0])
            for matrix in (normalized_candidate, normalized_reference)]
        ratios = eigvalsh(normalized_candidate, normalized_reference)
        report["variance_ratio_range"] = [float(ratios[0]), float(ratios[-1])]
        report["maximum_relative_mode_change"] = float(np.max(np.abs(ratios-1)))
        for name in names:
            actual, expected = candidate[name], reference[name]
            delta = actual-expected
            modes = eigvalsh(delta/normalization, normalized_reference)
            expected_diagonal = np.diag(expected)
            nonzero = expected_diagonal != 0
            diagonal_changes = np.diag(delta)[nonzero]/expected_diagonal[nonzero]
            report["components"][name] = dict(
                bitwise_equal=bool(np.array_equal(actual, expected)),
                relative_frobenius=float(np.linalg.norm(delta)
                                         /np.linalg.norm(expected)),
                max_change_over_total_diagonal=float(
                    np.max(np.abs(delta)/normalization)),
                max_relative_diagonal_change=float(
                    np.max(np.abs(diagonal_changes))),
                component_zero_diagonal_count=int(np.sum(~nonzero)),
                variance_change_range=[float(modes[0]), float(modes[-1])],
                maximum_total_variance_change=float(np.max(np.abs(modes))),
            )
    args.output.write_text(json.dumps(report, indent=2)+"\n")

    # A common reference-total normalization makes the component panels
    # comparable even when a component has tiny or zero off-diagonal entries.
    figure, axes = plt.subplots(2, 2, figsize=(11, 9), layout="constrained")
    probes = np.repeat(reference["rows"][:, 0], len(reference["coordinate"]))
    boundaries = np.flatnonzero(np.diff(probes))+1
    edges = np.concatenate(([0], boundaries, [len(diagonal)]))
    centers = (edges[:-1]+edges[1:])/2
    labels = [r"$\xi_+$", r"$\xi_-$", r"$\gamma_t$", r"$w$"]
    for axis, name in zip(axes.ravel(), names):
        delta = 100*(candidate[name]-reference[name])/normalization
        maximum = np.max(np.abs(delta))
        if maximum == 0:
            axis.set(xlim=(-0.5, len(diagonal)-0.5),
                     ylim=(len(diagonal)-0.5, -0.5))
            axis.text(0.5, 0.5, "Exactly zero difference", ha="center",
                      va="center", transform=axis.transAxes)
        else:
            picture = axis.imshow(delta, origin="upper", cmap="RdBu_r",
                                  vmin=-maximum, vmax=maximum,
                                  interpolation="nearest", rasterized=True)
            figure.colorbar(picture, ax=axis, shrink=0.85,
                            label="Difference / reference total rms product [%]")
        for boundary in boundaries:
            axis.axhline(boundary-0.5, color="0.3", linewidth=0.4)
            axis.axvline(boundary-0.5, color="0.3", linewidth=0.4)
        axis.set(xticks=centers, yticks=centers,
                 xticklabels=labels, yticklabels=labels, title=name.upper())
    minimum = first_record["mass_min"]
    reference_minimum = second_record["mass_min"]
    if minimum == reference_minimum:
        comparison = "Original versus expanded sigma-table domain"
        cutoffs = f"Both mass integrals start at {minimum:g} solar masses/h"
    else:
        comparison = "Common sigma-table domain"
        if first_record["interface_sha256"] != second_record["interface_sha256"]:
            comparison = "Original versus expanded sigma-table domain"
        cutoffs = (f"{minimum:g} minus {reference_minimum:g} solar masses/h "
                   "integration cutoff")
    if args.label is not None:
        comparison = args.label
    figure.suptitle(
        f"Complete LSST Y1 covariance · {comparison}\n{cutoffs}\n"
        "Every matrix entry retained; each panel has its own color scale",
        fontsize=13)
    figure.savefig(args.figure, dpi=180)
    figure.savefig(args.figure.with_suffix(".pdf"))
    print(json.dumps({"output": str(args.output),
                      "maximum_relative_mode_change": report[
                          "maximum_relative_mode_change"],
                      "both_positive_definite": True}))


if __name__ == "__main__":
    main()
