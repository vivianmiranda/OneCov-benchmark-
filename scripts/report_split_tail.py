"""Compare a cheap low-mass tail and extrapolated native partial sums.

Finite-k I11 is tested both before and after its existing completion.
The known zero-k bias integral is one; predicting that identity alone
does not demonstrate an improvement over the current halo prescription.
"""

import argparse
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

from common import sha256
from report_mass_cutoff import differences, read_run


def extrapolate(values):
    """Aitken on the last three sums; Wynn/Shanks on all supplied sums.

    Aitken removes a geometric remainder. Nearly equal sums provide no
    usable curvature, so report no estimate instead of dividing by noise.
    Extra arithmetic precision cannot recover accuracy lost in input sums.
    """
    curvature = values[-1]-2*values[-2]+values[-3]
    if abs(curvature) < 1e-12*max(abs(value) for value in values):
        return None, None
    aitken = values[-1]-(values[-1]-values[-2])**2/curvature
    with mp.workdps(40):
        table = mp.shanks([mp.mpf(value) for value in values], randomized=False)
        wynn = float(table[-1][-1])
    return float(aitken), wynn


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("single", "sequence", "reference", "refined", "output",
                 "figure"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.figure.exists():
        parser.error("choose new outputs")
    np.testing.assert_allclose(extrapolate([0.5, 0.75, 0.875]), [1, 1])
    data, records = {}, {}
    for name in ("single", "sequence", "reference", "refined"):
        data[name], records[name] = read_run(getattr(args, name))
    reference = data["reference"]
    if len({r["power_table_sha256"] for r in records.values()}) != 1:
        raise ValueError("power inputs differ")
    for run in data.values():
        for key in ("k", "first", "second", "redshift", "tail_k", "tail_power"):
            np.testing.assert_array_equal(run[key], reference[key])

    chosen = (reference["k"] <= 100,
              np.ones(len(reference["first"]), dtype=bool))
    checks = {}
    for name in ("single", "sequence", "refined"):
        checks[name+"_vs_reference"] = differences(
            data[name], reference, "m-20_n96", "m-20_n96", chosen)
    checks["cutoff4_vs20"] = differences(
        reference, reference, "m4_n96", "m-20_n96", chosen)

    exponents = [4, 0, -4, -8, -12, -16, -20]
    partials = {}
    for name in ("sequence", "refined"):
        rows = {row["key"]: row for row in records[name]["records"]}
        partials[name] = np.array([
            1-np.array(rows[f"m{exponent}_n96"]["missing_weight"])
            for exponent in exponents])

    # The final two partial sums are withheld from the geometric fit.
    # Its prediction of those actual integrals is distinct from the
    # extrapolated infinite-mass-range limit, whose bias sum is known.
    held_out = []
    for row in range(3):
        sums = partials["sequence"][:, row]
        delta = sums[4]-sums[3]
        ratio = delta/(sums[3]-sums[2])
        predicted = [sums[4]+delta*ratio,
                     sums[4]+delta*ratio+delta*ratio**2]
        held_out.append(dict(
            predicted=predicted,
            actual=partials["refined"][5:, row].tolist(),
            relative_error=(np.array(predicted)
                            /partials["refined"][5:, row]-1).tolist(),
        ))

    start = time.perf_counter()
    estimates = []
    finite_k = []
    for name in ("sequence", "refined"):
        run = data[name]
        for count in (3, 5, 7):
            bias = [extrapolate(partials[name][:count, row])
                    for row in range(3)]
            estimates.append(dict(
                run=name, partial_sum_count=count,
                minimum_mass=10.0**exponents[count-1],
                aitken=[value[0] for value in bias],
                wynn=[value[1] for value in bias],
            ))
            errors = []
            for row in range(3):
                for index in range(len(run["k"])):
                    # Remove each native minimum-profile completion to
                    # obtain actual resolved partial integrals of I11.
                    sums = [float(run[f"m{e}_n96_i11"][row, index]
                                  -run[f"m{e}_n96_correction"][row, index])
                            for e in exponents[:count]]
                    values = extrapolate(sums)
                    expected = reference["m-20_n96_i11"][row, index]
                    errors.append(np.array(values)/expected-1)
            finite_k.append(dict(
                run=name, partial_sum_count=count,
                maximum_fractional_error=np.max(np.abs(errors), axis=0).tolist(),
            ))
    extrapolation_seconds = time.perf_counter()-start

    # Once the existing completion is included, later low-mass extensions
    # can fall below the resolution of a curvature-based extrapolation.
    completed = np.array([data["sequence"][f"m{e}_n96_i11"]
                          for e in exponents])
    curvature = completed[-1]-2*completed[-2]+completed[-3]
    unresolved_curvature = np.abs(curvature) < 1e-12*np.max(
        np.abs(completed), axis=0)
    result = dict(
        redshift=reference["redshift"].tolist(), exponents=exponents,
        partial_bias_sums={k: v.tolist() for k, v in partials.items()},
        checks=checks, bias_limit_estimates=estimates,
        finite_k_resolved_I11_extrapolation=finite_k,
        held_out_partial_sums=held_out,
        completed_I11_unresolved_curvature_fraction=float(
            np.mean(unresolved_curvature)),
        extrapolation_seconds=extrapolation_seconds,
        mpmath_version=mp.__version__, runs=records,
        script_sha256=sha256(__file__),
    )
    args.output.write_text(json.dumps(result, indent=2)+"\n")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for row, redshift in enumerate(reference["redshift"]):
        label = f"z = {redshift:g}"
        axes[0].plot(exponents, 100*(1-partials["sequence"][:, row]),
                     "o-", label=label)
        selected = [item for item in estimates if item["run"] == "sequence"]
        limits = [np.log10(item["minimum_mass"]) for item in selected]
        axes[1].plot(limits, [100*(item["aitken"][row]-1)
                              for item in selected], "o-", label=label)
    axes[0].set(ylabel="Missing bias-weighted response [%]")
    axes[1].set(ylabel="Aitken estimate minus the known unit sum [%]")
    axes[1].axhline(0, color="black", linewidth=0.6)
    for axis in axes:
        axis.set_xlabel(r"Lowest included $\log_{10}(M/[M_\odot/h])$")
        axis.invert_xaxis()
        axis.grid(alpha=0.2)
    axes[0].legend()
    fig.suptitle("32 nodes per added four-decade interval; unchanged halo fits")
    fig.savefig(args.figure, dpi=180)
    fig.savefig(args.figure.with_suffix(".pdf"))
    print(f"Saved split-tail and partial-sum comparison: {args.output}")


if __name__ == "__main__":
    main()
