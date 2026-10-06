"""Summarize the controlled sigma and abundance diagnostics, with plots.

Run in the OneCov environment. Its native hmf filter evaluates sigma8
from the archived spectra; all other calculations are ratios of outputs
already saved by the two codes.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from hmf.density_field.filters import TopHat

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("onecov", type=Path)
    parser.add_argument("cocoa", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path, default=Path("results/figures"))
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file")
    data = np.load(args.onecov/"diagnostic.npz")
    cocoa = np.load(args.cocoa/"inputs.npz")
    report = json.loads((args.onecov/"report.json").read_text())
    cocoa_report = json.loads((args.cocoa/"report.json").read_text())
    for boost in (1, 2, 4):
        if not np.array_equal(data[f"cocoa_fft_boost{boost}"],
                              cocoa[f"sigma_boost{boost}"]):
            raise ValueError("the Cocoa export does not match these diagnostics")

    def maximum_ratio(left, right):
        return np.max(np.abs(left/right-1), axis=1).tolist()

    stages = ["onecov_native", "cocoa_power_native_onecov_radius",
              "cocoa_power_dense_onecov_radius", "cocoa_power_wide_onecov_radius",
              "shared_filter_stride1", "cocoa_fft_boost1"]
    transitions = {f"{a} -> {b}": maximum_ratio(data[b], data[a])
                   for a, b in zip(stages[:-1], stages[1:])}
    refinements = {f"{200*n} -> {400*n}": maximum_ratio(
        data[f"onecov_refine{n}"], data[f"onecov_refine{2*n}"])
        for n in (1, 2, 4, 8, 16, 32)}

    # Give the two native power spectra the same sigma8 only in this
    # explicit diagnostic. Do not change either production cosmology.
    sigma8_c = np.array([TopHat(cocoa["k_wide"], p).sigma(np.array([8.0]))[0]
                         for p in cocoa["power_wide"]])
    sigma8_o = data["onecov_sigma8_refine64"]
    common_amplitude = (data["cocoa_power_wide_onecov_radius"]
                        * (sigma8_o/sigma8_c)[:, None])
    z0 = {name: float(TopHat(cocoa["k_amplitude"], cocoa[name]).sigma(
        np.array([8.0]))[0]) for name in ("power_reader_z0", "power_camb_z0")}

    factors = data["abundance_multiplicity_at_same_peak"]
    residual = data["abundance_ratio"]/factors-1
    result = dict(
        redshifts=[row["z"] for row in report["rows"]],
        sigma_native_max_fractional=maximum_ratio(
            data["cocoa_fft_boost1"], data["onecov_native"]),
        sigma_controlled_steps=transitions,
        sigma_shared_max_fractional=maximum_ratio(
            data["cocoa_fft_boost1"], data["shared_filter_stride1"]),
        sigma_shared_boost4_max_fractional=maximum_ratio(
            data["cocoa_fft_boost4"], data["shared_filter_stride1"]),
        sigma_filter_refinement=maximum_ratio(
            data["shared_filter_stride2"], data["shared_filter_stride1"]),
        sigma_onecov_native_refinements=refinements,
        sigma_common_amplitude_max_fractional=maximum_ratio(
            common_amplitude, data["onecov_refine64"]),
        sigma8_cocoa_power=sigma8_c.tolist(), sigma8_onecov_fine=sigma8_o.tolist(),
        sigma8_z0=z0, abundance_factors=report["rows"],
        abundance_without_amplitude_range=[
            [float(row.min()), float(row.max())] for row in residual],
        source_reports=dict(onecov=report, cocoa=cocoa_report),
        inputs_sha256=dict(onecov=sha256(args.onecov/"diagnostic.npz"),
                           cocoa=sha256(args.cocoa/"inputs.npz")),
        script_sha256=sha256(__file__))
    args.output.write_text(json.dumps(result, indent=2)+"\n")

    figure, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    mass = data["mass"]
    for row, z in enumerate(result["redshifts"]):
        color = f"C{row}"
        axes[0, 0].semilogx(mass, 100*(data["cocoa_fft_boost1"][row]
                                     /data["onecov_native"][row]-1),
                            color=color, label=f"z={z:g}")
        axes[0, 1].semilogx(mass, 100*(data["cocoa_fft_boost1"][row]
                                     /data["shared_filter_stride1"][row]-1),
                            color=color)
        axes[1, 0].semilogx(mass, 100*(data["abundance_ratio"][row]-1), color=color)
        axes[1, 0].axhline(100*(factors[row, 0]-1), color=color, linestyle="--")
        axes[1, 1].semilogx(mass, 100*residual[row], color=color)
    titles = ("Native mass rms fluctuation", "Same power and radius: sigma",
              "Native abundance; dashed = normalization only",
              "Abundance after removing normalization factor")
    for axis, title in zip(axes.ravel(), titles):
        axis.axhline(0, color="0.6", linewidth=0.7)
        axis.set(title=title, xlabel=r"$M\ [M_\odot/h]$", ylabel="Difference [%]")
    axes[0, 0].legend()
    figure.suptitle("CoCoA / OneCov − 1 · controlled halo-ingredient comparisons")
    args.figures.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(args.figures/f"sigma_abundance.{extension}", dpi=180)
    plt.close(figure)
    print(json.dumps({k: result[k] for k in ("sigma8_z0",
        "sigma_native_max_fractional", "sigma_shared_max_fractional",
        "sigma_common_amplitude_max_fractional",
        "abundance_without_amplitude_range")}, indent=2))


if __name__ == "__main__":
    main()
