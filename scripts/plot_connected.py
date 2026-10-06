"""Plot native halo terms and shared-input cNG projection; never timings."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trispectrum", type=Path)
    parser.add_argument("projection", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = np.load(args.trispectrum/"comparison.npz")
    diagonal = data["first"] == data["second"]
    colors = ["#2465a8", "#bf541b", "#25814d"]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    for term, ax in enumerate(axes.flat[:4]):
        for row, (z, color) in enumerate(zip(data["redshift"], colors)):
            ax.loglog(data["k"], data["native"][row, term, diagonal],
                      color=color, label=f"z={z:g}")
            ax.loglog(data["k"], data["onecov"][row, term, diagonal],
                      color=color, linestyle="--")
        ax.set(title=f"{term+1}h: native diagonal", xlabel=r"$k$ [h/Mpc]",
               ylabel=r"$T(k,k)$ [(Mpc/h)$^9$]")
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=9)
    size = len(data["k"])
    residual = np.zeros((size, size))
    residual[data["first"], data["second"]] = 100*(
        data["shared"][1, 1]/data["onecov"][1, 1]-1)
    residual[data["second"], data["first"]] = residual[
        data["first"], data["second"]]
    im = axes[1, 1].imshow(residual, origin="lower", cmap="magma", vmin=0)
    fig.colorbar(im, ax=axes[1, 1], label="CoCoA / OneCov − 1 [%]")
    axes[1, 1].set(title="2h with shared inputs, z=0.5",
                   xlabel="k index (0.001 → 10 h/Mpc)", ylabel="q index")
    for row, (z, color) in enumerate(zip(data["redshift"], colors)):
        ratio = data["native"][row].sum(axis=0)/data["onecov"][row].sum(axis=0)
        axes[1, 2].semilogx(data["k"], ratio[diagonal], color=color,
                           label=f"z={z:g}")
    axes[1, 2].axhline(1, color="black", linewidth=0.7)
    axes[1, 2].set(title="Native sum: different prescriptions",
                   xlabel=r"$k$ [h/Mpc]", ylabel="CoCoA / OneCov")
    axes[1, 2].grid(alpha=0.2)
    fig.suptitle("Matter trispectrum · CoCoA solid, OneCov dashed\n"
                 "Native models retain their halo fits, damping and corner treatments")
    for suffix in ("png", "pdf"):
        fig.savefig(args.output/f"trispectrum_terms.{suffix}", dpi=160)
    plt.close(fig)

    data = np.load(args.projection/"comparison.npz")
    reference = data["onecov"]
    rms = np.sqrt(reference.diagonal())
    residual = (data["cocoa"]-reference)/rms[:, None]/rms[None, :]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    axes[0].loglog(data["ell"], data["cocoa"].diagonal(), label="CoCoA")
    axes[0].loglog(data["ell"], reference.diagonal(), "o", fillstyle="none",
                   label="OneCovariance")
    axes[0].set(xlabel=r"$\ell$", ylabel=r"$C^{\rm cNG}_{\ell\ell}$",
                title="Shared matter tables and radial weights")
    axes[0].legend()
    limit = max(np.max(np.abs(residual))*1e15, 1)
    im = axes[1].imshow(residual*1e15, origin="lower", cmap="RdBu_r",
                        vmin=-limit, vmax=limit)
    fig.colorbar(im, ax=axes[1], label=r"$10^{15}\Delta C_{ij}/\sqrt{C_{ii}C_{jj}}$")
    axes[1].set(xlabel="Multipole index", ylabel="Multipole index",
                title="All 64 entries, including crossed multipoles")
    fig.suptitle("LSST Y1 source bin 3 · Connected non-Gaussian projection")
    for suffix in ("png", "pdf"):
        fig.savefig(args.output/f"connected_projection.{suffix}", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
