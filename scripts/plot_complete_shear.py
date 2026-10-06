"""Plot CoCoA minus OneCovariance for every entry of a complete shear pilot.

All four panels use the same OneCov total-variance normalization. Separate
colour ranges reveal small component differences without dividing by tiny
off-diagonal entries. Inputs must come from complete_shear.py on identical
source distributions, cosmology, density and measured multipoles or annuli.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import cholesky, eigvalsh
from threadpoolctl import threadpool_limits

from common import sha256


def read_result(folder):
    """Read four observable-by-observable covariances and check their hash."""
    report = json.loads((folder/"report.json").read_text())
    filename = folder/"covariance.npz"
    if sha256(filename) != report["covariance_sha256"]:
        raise ValueError(f"matrix archive changed: {filename}")
    return dict(np.load(filename, allow_pickle=False)), report


def main():
    """Check every component and total, then save PNG, PDF and diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cocoa", type=Path)
    parser.add_argument("onecov", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.figure.exists():
        parser.error("choose new output and figure paths")
    cocoa, cocoa_record = read_result(folder=args.cocoa)
    onecov, onecov_record = read_result(folder=args.onecov)
    for key in ("input_manifest_sha256", "source_bin_1based", "delta_ell"):
        if cocoa_record[key] != onecov_record[key]:
            raise ValueError(f"comparison inputs differ: {key}")
    space = cocoa_record.get("space", "fourier")
    if space != onecov_record.get("space", "fourier"):
        raise ValueError("both calculations must use the same observable space")
    coordinate = "theta_edges_arcmin" if space == "real" else "ell"
    np.testing.assert_allclose(cocoa[coordinate], onecov[coordinate], rtol=1.e-13)
    names = ("gaussian", "ssc", "cng", "total")
    diagonal = np.diag(onecov["total"])
    if np.any(diagonal <= 0):
        raise ValueError("OneCov total has a nonpositive variance")
    rms_product = np.sqrt(diagonal[:, None]*diagonal[None, :])

    # A common total scale expresses each component difference in units
    # of the measured uncertainty, including zero off-diagonal entries.
    # Positivity applies to the total, not separately to every halo term.
    for result in (cocoa, onecov):
        for name in names:
            matrix = result[name]
            if matrix.shape != rms_product.shape or not np.isfinite(matrix).all():
                raise ValueError(f"invalid {name} component")
            asymmetry = np.max(np.abs(matrix-matrix.T)/rms_product)
            if asymmetry > 1.e-4:
                raise ValueError(f"asymmetric {name} component")
        summed = result["gaussian"]+result["ssc"]+result["cng"]
        np.testing.assert_array_equal(result["total"], summed)

    report = {
        "scope": cocoa_record["scope"],
        "difference": "CoCoA minus OneCovariance",
        "normalization": "OneCov total diagonal rms product",
        "ndata": len(diagonal),
        "masked_entries": 0,
        "cocoa": cocoa_record,
        "onecov": onecov_record,
        "script_sha256": sha256(__file__),
        "components": {},
    }
    # Native double projections may leave small antisymmetric residuals.
    # Use symmetric parts only for variance eigenvalues, never for the
    # plotted differences or archived matrices. Report the residuals too.
    reference = (onecov["total"]+onecov["total"].T)/(2*rms_product)
    candidate = (cocoa["total"]+cocoa["total"].T)/(2*rms_product)
    report["eigenvalue_convention"] = "symmetric parts of native matrices"
    with threadpool_limits(limits=1, user_api="blas"):
        cholesky(a=reference, lower=True)
        cholesky(a=candidate, lower=True)
        report["both_positive_definite"] = True
        ratios = eigvalsh(a=candidate, b=reference)
        report["total_variance_ratio_range"] = [float(ratios[0]), float(ratios[-1])]
        for name in names:
            delta = cocoa[name]-onecov[name]
            modes = eigvalsh(a=(delta+delta.T)/(2*rms_product), b=reference)
            report["components"][name] = {
                "max_total_variance_scaled_difference_percent":
                    float(100*np.max(np.abs(delta)/rms_product)),
                "max_relative_diagonal_difference_percent":
                    float(100*np.max(np.abs(np.diag(delta)/np.diag(onecov[name])))),
                "total_variance_mode_change_range_percent":
                    [float(100*modes[0]), float(100*modes[-1])],
                "max_asymmetry_total_rms_scaled": {
                    "cocoa": float(np.max(np.abs(cocoa[name]-cocoa[name].T)/rms_product)),
                    "onecov": float(np.max(np.abs(onecov[name]-onecov[name].T)/rms_product)),
                },
            }
    args.output.write_text(json.dumps(report, indent=2)+"\n")

    plt.rcParams.update({
        "font.family": "STIXGeneral",
        "mathtext.fontset": "stix",
        "axes.labelsize": 24,
        "axes.titlesize": 24,
        "xtick.labelsize": 19,
        "ytick.labelsize": 19,
        "axes.linewidth": 1.0,
    })
    figure, axes = plt.subplots(2, 2, figsize=(14, 12), layout="constrained")
    ticks = np.array([0, 24, 49, 74, 99])
    labels = None
    if space == "fourier":
        labels = onecov["ell"][ticks].astype(int)
    else:
        ticks = [3.5, 11.5]
        labels = [r"$\xi_+$", r"$\xi_-$"]
    titles = ("Gaussian", "SSC", "Connected non-Gaussian", "Total")
    for axis, name, title in zip(axes.ravel(), names, titles):
        delta_percent = 100*(cocoa[name]-onecov[name])/rms_product
        maximum = float(np.max(np.abs(delta_percent)))
        picture = axis.imshow(
            delta_percent, origin="upper", cmap="RdBu_r",
            vmin=-maximum, vmax=maximum, interpolation="nearest",
            rasterized=True,
        )
        colourbar = figure.colorbar(picture, ax=axis, shrink=0.88, pad=0.03)
        colourbar.set_label("Difference / OneCov total rms product [%]", size=18)
        colourbar.ax.tick_params(labelsize=16)
        axis.set(
            title=title,
            xlabel=r"$\ell$" if space == "fourier" else "Angular bins",
            ylabel=r"$\ell'$" if space == "fourier" else "Angular bins",
            xticks=ticks, yticks=ticks, xticklabels=labels, yticklabels=labels,
        )
        if space == "real":
            axis.axhline(7.5, color="0.35", lw=0.8)
            axis.axvline(7.5, color="0.35", lw=0.8)
    description = "100 × 100 Fourier covariance · native halo prescriptions"
    if space == "real":
        description = ("Real-space covariance · 8 bins each: 2.5–250 arcmin\n"
                       "CoCoA full sky · OneCovariance flat sky")
    figure.suptitle(
        "LSST Y1 source bin 3 · CoCoA − OneCovariance\n"
        +description+"\n"
        "Every matrix entry retained; each panel has its own colour scale",
        fontsize=23,
    )
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")
    figure.savefig(args.figure.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(report["components"], indent=2))
    print("Both total symmetric parts are positive definite.")


if __name__ == "__main__":
    main()
