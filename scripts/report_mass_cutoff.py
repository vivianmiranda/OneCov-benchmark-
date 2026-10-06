"""Summarize the unchanged-model cutoff experiment and plot its physics.

All ratios use saved production-kernel outputs. Times cover halo moments
only; no timing plot or full-covariance speed estimate is produced.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import sha256


def read_run(folder):
    record = json.loads((folder/"report.json").read_text())
    if sha256(folder/"ingredients.npz") != record["ingredients_sha256"]:
        raise ValueError(f"ingredients changed: {folder}")
    if sha256(folder/"power_tables.npz") != record["power_table_sha256"]:
        raise ValueError(f"power changed: {folder}")
    data = dict(np.load(folder/"ingredients.npz"))
    for name, values in data.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"nonfinite {name} in {folder}")
    return data, record


def differences(left, right, prefix_left, prefix_right, chosen):
    """Largest absolute fractional changes over the chosen k or pair set."""
    result = {}
    for name in ("i11", "moments", "response", "terms"):
        mask = chosen[1] if name in ("moments", "terms") else chosen[0]
        lhs = left[prefix_left+"_"+name][..., mask]
        rhs = right[prefix_right+"_"+name][..., mask]
        if np.any(rhs == 0):
            raise ValueError(f"zero denominator in {name}")
        result[name] = float(np.max(np.abs(lhs/rhs-1)))
        if name == "terms":
            result["each_term"] = np.max(
                np.abs(lhs/rhs-1), axis=(0, 2)).tolist()
            result["total_trispectrum"] = float(np.max(np.abs(
                lhs.sum(axis=1)/rhs.sum(axis=1)-1)))
        if name == "moments":
            result["each_moment"] = np.max(
                np.abs(lhs/rhs-1), axis=(1, 2)).tolist()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--power-check", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file")
    fine, fine_record = read_run(args.fine)
    coarse, coarse_record = read_run(args.coarse)
    native, native_record = read_run(args.native)
    fingerprints = [r["power_table_sha256"] for r in
                    (fine_record, coarse_record, native_record)]
    if len(set(fingerprints)) != 1:
        raise ValueError("runs used different supplied power tables")
    power_check = json.loads(args.power_check.read_text())
    if power_check["source_ingredients_sha256"] != fine_record[
            "ingredients_sha256"]:
        raise ValueError("power check does not refer to the fine run")

    k = fine["k"]
    first, second = fine["first"], fine["second"]
    comparisons = {}
    for maximum in (10, 100):
        points = k <= maximum
        chosen = (points, points[first] & points[second])
        transitions = {}
        for low, high in ((6, 4), (4, 2), (6, 2)):
            transitions[f"1e{low}_vs_1e{high}"] = differences(
                fine, fine, f"m{low}_n512", f"m{high}_n512", chosen)
        checks = {}
        for exponent in (6, 4, 2):
            checks[f"quadrature_m{exponent}"] = differences(
                fine, fine, f"m{exponent}_n256", f"m{exponent}_n512", chosen)
            checks[f"tables_m{exponent}"] = differences(
                coarse, fine, f"m{exponent}_n256", f"m{exponent}_n256", chosen)
        checks["native_vs_wide_table"] = differences(
            native, fine, "m6_n256", "m6_n256", chosen)
        # Common table errors cancel strongly in a cutoff contrast. Check
        # that cancellation explicitly rather than comparing error scales.
        contrast_checks = {}
        for low, high in ((6, 4), (4, 2), (6, 2)):
            changes = {}
            for name in ("i11", "moments", "response", "terms"):
                mask = chosen[1] if name in ("moments", "terms") else chosen[0]
                old = (coarse[f"m{low}_n256_{name}"][..., mask]
                       /coarse[f"m{high}_n256_{name}"][..., mask]-1)
                new = (fine[f"m{low}_n512_{name}"][..., mask]
                       /fine[f"m{high}_n512_{name}"][..., mask]-1)
                changes[name] = float(np.max(np.abs(old-new)))
            contrast_checks[f"1e{low}_vs_1e{high}"] = changes
        comparisons[str(maximum)] = dict(
            cutoff_changes=transitions, numerical_checks=checks,
            cutoff_contrast_stability=contrast_checks)

    timing = []
    for item in fine_record["records"]:
        if item["nquad"] != 512:
            continue
        samples = np.asarray(item["seconds"])
        timing.append(dict(mass_min=item["mass_min"],
                           mass_panels=item["mass_panels"],
                           mean_seconds=float(samples.mean()),
                           std_seconds=float(samples.std(ddof=1))))
    record = dict(
        mass_cutoffs=[1e6, 1e4, 1e2], redshift=fine["redshift"].tolist(),
        missing_weights=[fine[f"m{m}_n512_missing"].tolist() for m in (6, 4, 2)],
        comparisons=comparisons, timing=timing, power_coverage=power_check,
        moment_order=["I02", "I12", "I13(K,Q,Q)", "I13(K,K,Q)", "I04"],
        term_order=["1h", "2h13", "2h22", "3h", "4h"],
        scope="CoCoA cutoff diagnostic; same model, no survey projection",
        source_reports=dict(fine=fine_record, coarse=coarse_record,
                            native=native_record),
        script_sha256=sha256(__file__),
    )
    args.output.write_text(json.dumps(record, indent=2)+"\n")

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for row, z in enumerate(fine["redshift"]):
        missing = [fine[f"m{m}_n512_missing"][row] for m in (2, 4, 6)]
        axes[0].semilogx([1e2, 1e4, 1e6], 100*np.array(missing),
                        marker="o", color=f"C{row}", label=f"z={z:g}")
        for exponent, style in ((6, "-"), (4, "--")):
            ratio = (fine[f"m{exponent}_n512_i11"][row]
                     /fine["m2_n512_i11"][row]-1)
            axes[1].plot(k, 100*np.abs(ratio), color=f"C{row}",
                         linestyle=style)
    axes[0].set(xlabel=r"Minimum halo mass $[M_\odot/h]$",
                ylabel="Missing bias-weighted response [%]",
                title="Correction weight decreases slowly")
    axes[0].legend()
    axes[1].set(xscale="log", yscale="log", ylim=(1e-8, 1),
                xlabel=r"$k\ [h/\mathrm{Mpc}]$",
                ylabel=r"Corrected $I_{11}$: difference from lowest cutoff [%]",
                title=r"Solid: $10^6$; dashed: $10^4$ vs $10^2\ M_\odot/h$")
    for axis in axes:
        axis.grid(alpha=0.2)
    figure.suptitle("CoCoA · same halo model and continued power spectrum")
    destination = args.output.parent/"figures"
    destination.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(destination/f"mass_cutoff.{extension}", dpi=180)
    print(f"Saved cutoff report: {args.output}")


if __name__ == "__main__":
    main()
