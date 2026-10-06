"""Show the shared-input SSC projection and its numerical refinements."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison", type=Path)
    parser.add_argument("--refinements", type=Path, nargs="*", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = np.load(args.comparison / "matrices.npz")
    ell = data["ell"]
    onecov = data["onecov"]
    cocoa = data["cocoa"]
    rms = np.sqrt(np.diag(onecov))
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")

    axes[0, 0].loglog(ell, np.diag(cocoa), label="CoCoA: shared inputs")
    axes[0, 0].loglog(ell[::5], np.diag(onecov)[::5], "o", fillstyle="none",
                       label="OneCovariance")
    axes[0, 0].set(xlabel=r"$\ell$", ylabel=r"$\mathrm{Cov}_{\rm SSC}(C_\ell,C_\ell)$",
                   title="Shear SSC diagonal")
    axes[0, 0].legend()

    residual = (cocoa-onecov)/rms[:, None]/rms[None, :]
    image = axes[0, 1].imshow(residual*1.e15, origin="lower", cmap="RdBu_r")
    axes[0, 1].set(title="CoCoA − OneCov: all matrix entries",
                   xlabel="Multipole sample", ylabel="Multipole sample")
    figure.colorbar(image, ax=axes[0, 1], label=r"$10^{15}\Delta C_{ij}/\sqrt{C_{ii}C_{jj}}$")

    # Only the footprint changes in this panel. Shared response and power
    # arrays keep the comparison separate from native halo-model choices.
    cap_change = np.diag(data["analytic_cap"])/np.diag(onecov)-1
    axes[1, 0].semilogx(ell, cap_change*100)
    axes[1, 0].set(title="Analytic versus pixelized spherical cap",
                   xlabel=r"$\ell$", ylabel="SSC variance change [%]")

    for folder in args.refinements:
        refined = np.load(folder / "matrices.npz")
        if not np.array_equal(refined["ell"], ell):
            raise ValueError("convergence plots require identical multipoles")
        change = np.diag(refined["onecov"])/np.diag(onecov)-1
        settings = json.loads((folder/"comparison.json").read_text())["export"]
        label = (f"{settings['radial_nodes']} radial nodes; "
                 f"{settings['response_redshifts']} response rows")
        axes[1, 1].semilogx(ell, change*100, label=label)
    axes[1, 1].set(title="OneCov refinement versus the base grid",
                   xlabel=r"$\ell$", ylabel="SSC variance change [%]")
    if args.refinements:
        axes[1, 1].legend(fontsize=8)
    else:
        correlation = onecov/rms[:, None]/rms[None, :]
        image = axes[1, 1].imshow(correlation, origin="lower", vmin=-1,
                                   vmax=1, cmap="RdBu_r")
        axes[1, 1].set(title="OneCov SSC correlations",
                       xlabel="Multipole sample", ylabel="Multipole sample")
        figure.colorbar(image, ax=axes[1, 1], label="Correlation")
    figure.suptitle("LSST Y1 source bin 3 · SSC projection with shared ingredients")
    args.output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(args.output/f"ssc_projection.{extension}", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
