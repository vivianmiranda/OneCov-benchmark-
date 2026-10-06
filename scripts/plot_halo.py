"""Plot native CoCoA/OneCov halo differences, retaining each prescription."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("onecov", type=Path)
    parser.add_argument("cocoa", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    native = np.load(args.onecov/"inputs.npz")
    cocoa = np.load(args.cocoa/"cocoa.npz")
    quantities = (
        ("sigma", r"Mass rms fluctuation $\sigma(M)$", True),
        ("dndlnm", r"Abundance $dn/d\ln M$", True),
        ("bias", "Native halo bias", True),
        ("concentration", "Native concentration", True),
        ("linear", r"Linear power $P_{\rm lin}$", False),
        ("i11", r"Moment $I_1^1$", False),
        ("i02", r"Moment $I_2^0(k,k)$", False),
        ("i12", r"Moment $I_2^1(k,k)$", False),
        ("response", r"Native response $dP/d\delta_b$", False))
    figure, axes = plt.subplots(3, 3, figsize=(12, 10), layout="constrained")
    for axis, (name, title, mass_axis) in zip(axes.ravel(), quantities):
        for row, redshift in enumerate(native["redshift"]):
            x = native["mass" if mass_axis else "k"][row]
            selected = ((x >= 1.e10) & (x <= 1.e15) if mass_axis
                        else (x >= 0.001) & (x <= 10))
            difference = cocoa[name][row, selected]/native[name][row, selected]-1
            axis.semilogx(x[selected], difference*100, label=f"z={redshift:g}")
        axis.axhline(0, color="0.5", linewidth=0.7)
        axis.set(title=title, ylabel="CoCoA / OneCov − 1 [%]",
                 xlabel=r"$M\ [M_\odot/h]$" if mass_axis else r"$k\ [h/{\rm Mpc}]$")
    axes[0, 0].legend()
    figure.suptitle("Native halo prescriptions · matched LSST Y1 cosmology")
    args.output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(args.output/f"halo_ingredients.{extension}", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
