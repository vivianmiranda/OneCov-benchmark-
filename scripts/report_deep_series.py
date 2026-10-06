"""Check deeper low-mass partial sums against their known unit bias limit.

Use equal four-decade increments through 1e-48. The final 1e-50 integral
is a held-out check; its two-decade increment is not fed to an algorithm
that assumes equally spaced members of the sequence.
"""

import argparse
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import sha256
from report_mass_cutoff import differences, read_run
from report_split_tail import extrapolate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("coarse", "fine", "limits", "old", "output", "figure"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.figure.exists():
        parser.error("choose new output files")
    data, records = {}, {}
    for name in ("coarse", "fine"):
        data[name], records[name] = read_run(getattr(args, name))
    limits, limit_record = read_run(args.limits)
    if limit_record["boost"] != records["fine"]["boost"]:
        raise ValueError("full-fit normalization needs the fine table settings")
    if limit_record["diagnostic_build"] != records["fine"]["diagnostic_build"]:
        raise ValueError("normalization and partial sums use different builds")
    mass_limit = limits["full_mass_limit"]
    old, old_record = read_run(args.old)
    if old_record["power_table_sha256"] != records["fine"]["power_table_sha256"]:
        raise ValueError("old-cutoff control used different power inputs")
    for name in ("k", "redshift"):
        np.testing.assert_array_equal(old[name], data["fine"][name])
    old_row = next(row for row in old_record["records"]
                   if row["key"] == "m6_n256")
    if records["coarse"]["power_table_sha256"] != records["fine"][
            "power_table_sha256"]:
        raise ValueError("power inputs differ")
    for name in ("k", "first", "second", "redshift", "tail_k", "tail_power"):
        np.testing.assert_array_equal(data["coarse"][name], data["fine"][name])
    exponents = list(range(4, -49, -4))
    partials, estimates, finite_k, held_out = {}, [], [], []
    mass_partials, mass_estimates, deepest_mass = {}, [], {}
    start = time.perf_counter()
    for name in ("coarse", "fine"):
        rows = {r["key"]: r for r in records[name]["records"]}
        partials[name] = np.array([
            1-np.array(rows[f"m{e}_n96"]["missing_weight"])
            for e in exponents])
        mass_partials[name] = np.array([
            rows[f"m{e}_n96"]["resolved_mass"] for e in exponents])
        deepest_mass[name] = rows["m-50_n96"]["resolved_mass"]
        for count in range(3, len(exponents)+1):
            # Wynn's table uses an odd number of inputs for its final
            # limit estimate. Drop the oldest sum when necessary, so
            # every estimate includes the newest, deeper integral.
            window = min(count, 7)
            if window % 2 == 0:
                window -= 1
            whole = count if count % 2 else count-1
            bias = []
            mass = []
            for row in range(3):
                sums = partials[name][:count, row]
                aitken, wynn = extrapolate(sums[-window:])
                _, higher_order = extrapolate(sums[-whole:])
                bias.append([aitken, wynn, higher_order])
                # Repeat on ordinary mass weights, without halo bias.
                # Its limit belongs to the retained fitted model; it is
                # not forced to one by a numerical extrapolation.
                sums = mass_partials[name][:count, row]
                aitken, wynn = extrapolate(sums[-window:])
                _, higher_order = extrapolate(sums[-whole:])
                mass.append([aitken, wynn, higher_order])
            estimates.append(dict(
                run=name, exponent=exponents[count-1],
                wynn_window=window, wynn_full_window=whole,
                estimates=np.array(bias).tolist(),
            ))
            mass_estimates.append(dict(
                run=name, exponent=exponents[count-1],
                estimates=np.array(mass).tolist(),
                fractional_error=(np.array(mass)/mass_limit[:, None]-1).tolist(),
            ))

        # The geometric remainder fitted at -40,-44,-48 predicts the
        # additional half interval to -50. This point was not fitted.
        sums = partials[name][-3:]
        first = sums[1]-sums[0]
        last = sums[2]-sums[1]
        ratio = last/first
        predicted = sums[2]+last*ratio/(1-ratio)*(1-np.sqrt(ratio))
        actual = 1-np.array(rows["m-50_n96"]["missing_weight"])
        held_out.append(dict(run=name, predicted=predicted.tolist(),
                             actual=actual.tolist(),
                             relative_error=(predicted/actual-1).tolist()))

        # Remove the already-applied completion before extrapolating.
        # Otherwise the algorithm mainly sees the cancellation that the
        # current prescription has already supplied analytically.
        for count in range(3, len(exponents)+1):
            whole = count if count % 2 else count-1
            window = min(whole, 7)
            errors = []
            for row in range(3):
                for index in range(len(data[name]["k"])):
                    sums = [float(data[name][f"m{e}_n96_i11"][row, index]
                                  -data[name][f"m{e}_n96_correction"][row, index])
                            for e in exponents[:count]]
                    aitken, wynn = extrapolate(sums[-window:])
                    _, higher_order = extrapolate(sums[-whole:])
                    expected = data["fine"]["m-50_n96_i11"][row, index]
                    errors.append(np.array([aitken, wynn, higher_order])
                                  /expected-1)
            finite_k.append(dict(
                run=name, exponent=exponents[count-1],
                maximum_fractional_error=np.max(np.abs(errors), axis=0).tolist(),
            ))
    seconds = time.perf_counter()-start
    chosen = (data["fine"]["k"] <= 100,
              np.ones(len(data["fine"]["first"]), dtype=bool))
    result = dict(
        redshift=data["fine"]["redshift"].tolist(), exponents=exponents,
        estimate_order=["aitken", "wynn_last7", "wynn_highest_available_order"],
        partial_bias_sums={key: value.tolist() for key, value in partials.items()},
        partial_mass_sums={key: value.tolist()
                           for key, value in mass_partials.items()},
        mass_at_cutoff50=deepest_mass, full_mass_limit=mass_limit.tolist(),
        mass_limit_estimates=mass_estimates, normalization_run=limit_record,
        old_cutoff6=dict(
            run=old_record, missing_weight=old_row["missing_weight"],
            maximum_I11_fractional_error=float(np.max(np.abs(
                old["m6_n256_i11"]/data["fine"]["m-50_n96_i11"]-1))),
        ),
        bias_limit_estimates=estimates, finite_k_extrapolation=finite_k,
        held_out_mass50=held_out, extrapolation_seconds=seconds,
        corrected_cutoff20_vs50=differences(
            data["fine"], data["fine"], "m-20_n96", "m-50_n96", chosen),
        refinement=differences(
            data["coarse"], data["fine"], "m-50_n96", "m-50_n96", chosen),
        runs=records, script_sha256=sha256(__file__),
    )
    args.output.write_text(json.dumps(result, indent=2)+"\n")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for row, redshift in enumerate(result["redshift"]):
        axes[0].semilogy(exponents, 100*(1-partials["coarse"][:, row]),
                         "o-", label=f"z = {redshift:g}")
    selected = [row for row in estimates if row["run"] == "coarse"]
    for index, method in enumerate(("Aitken", "Wynn, last 7", "Wynn, full")):
        error = [100*abs(row["estimates"][2][index]-1) for row in selected]
        axes[1].semilogy([row["exponent"] for row in selected], error,
                         "o-", label=method)
    axes[0].set_ylabel("Missing bias-weighted response [%]")
    axes[1].set_ylabel("Limit-estimate error at z = 1 [%]")
    for axis in axes:
        axis.set_xlabel(r"Lowest included $\log_{10}(M/[M_\odot/h])$")
        axis.invert_xaxis()
        axis.grid(alpha=0.2)
        axis.legend()
    fig.suptitle("Deeper partial sums: unchanged extrapolated halo model")
    fig.savefig(args.figure, dpi=180)
    fig.savefig(args.figure.with_suffix(".pdf"))
    print(f"Saved deeper series comparison: {args.output}")


if __name__ == "__main__":
    main()
