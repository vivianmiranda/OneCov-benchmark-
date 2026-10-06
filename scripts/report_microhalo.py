"""Summarize the isolated 1e-3 Msun/h cutoff and its numerical checks.

Read production-kernel outputs at three sigma-integral endpoints and two
internal table resolutions. Keep integration convergence separate from the
physical uncertainty in the unchanged log-log power extrapolation.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import sha256
from report_mass_cutoff import differences, read_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("k5", "k7", "k9", "refined", "power_check", "output"):
        parser.add_argument("--"+name.replace("_", "-"), type=Path,
                            required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.figure.exists():
        parser.error("choose new output and figure paths")

    runs = {}
    provenance = {}
    for name in ("k5", "k7", "k9", "refined"):
        runs[name], provenance[name] = read_run(getattr(args, name))
    fine = runs["refined"]
    for data in runs.values():
        for name in ("k", "mass", "redshift", "tail_k", "tail_power"):
            np.testing.assert_array_equal(data[name], fine[name])
    fingerprints = [value["power_table_sha256"] for value in provenance.values()]
    if len(set(fingerprints)) != 1:
        raise ValueError("supplied CAMB tables differ")
    power = json.loads(args.power_check.read_text())
    if power["source_ingredients_sha256"] != provenance["refined"][
            "ingredients_sha256"]:
        raise ValueError("the power check must use the refined run")

    # Each integration cutoff keeps exactly the same retained mass panels.
    # F and B refer to ordinary mass and mass-weighted bias, respectively.
    exponents = (4, 2, 0, -3)
    rows = []
    for exponent in exponents:
        key = f"m{exponent}_n512"
        record = next(value for value in provenance["refined"]["records"]
                      if value["key"] == key)
        index = np.argmin(np.abs(fine["mass"]/(10.0**exponent)-1))
        np.testing.assert_allclose(fine["mass"][index], 10.0**exponent,
                                   rtol=1e-13)
        radius = (3*10.0**exponent/(4*np.pi*float(fine["rho"])))**(1/3)
        rows.append(dict(
            mass_min=10.0**exponent,
            mass_panels=record["mass_panels"],
            radius_mpc_h=radius,
            inverse_radius_hmpc=1/radius,
            missing_response=record["missing_weight"],
            resolved_mass=record["resolved_mass"],
            minimum_peak_height=(1.686/fine["sigma"][:, index]).tolist(),
        ))

    # Compare all k pairs, including unequal pairs, in two physical ranges.
    # Report every quadrature level rather than skipping the 128-node rule.
    checks = {}
    for maximum in (10, 100):
        points = fine["k"] <= maximum
        pairs = points[fine["first"]] & points[fine["second"]]
        chosen = (points, pairs)
        quadrature = {}
        for exponent in exponents:
            for nodes in (96, 128, 256):
                key = f"m{exponent}_n{nodes}"
                quadrature[key] = differences(
                    runs["k9"], runs["k9"], key,
                    f"m{exponent}_n512", chosen)
        tables = {}
        for exponent in exponents:
            key = f"m{exponent}_n512"
            tables[key] = differences(runs["k9"], fine, key, key, chosen)
        cutoffs = {}
        for exponent in (2, 0, -3):
            cutoffs[str(exponent)] = differences(
                fine, fine, "m4_n512", f"m{exponent}_n512", chosen)
        # Absolute table errors can exceed a small cutoff effect. Check
        # the effect itself in both table resolutions, with a common grid
        # for each ratio, before claiming that the cutoff change is resolved.
        contrast = {}
        for name in ("i11", "moments", "response", "terms"):
            selected = pairs if name in ("moments", "terms") else points
            ratios = []
            for data in (runs["k9"], fine):
                ratios.append(data["m4_n512_"+name][..., selected]
                              /data["m-3_n512_"+name][..., selected]-1)
            contrast[name] = float(np.max(np.abs(ratios[0]-ratios[1])))
        checks[str(maximum)] = dict(
            mass_quadrature=quadrature,
            table_boost4_vs8=tables,
            mass_cutoff4_vs_lower=cutoffs,
            cutoff_contrast_boost4_vs8=contrast,
        )

    sigma_checks = {}
    for left, right in (("k5", "k7"), ("k7", "k9"), ("k9", "refined")):
        difference = runs[left]["sigma"]/runs[right]["sigma"]-1
        sigma_checks[left+"_vs_"+right] = dict(
            maximum_fractional=float(np.max(np.abs(difference))),
            at_minimum_mass=difference[:, 0].tolist(),
        )
    report = dict(
        redshift=fine["redshift"].tolist(),
        cutoff_rows=rows,
        checks=checks,
        sigma_checks=sigma_checks,
        shared_power_filter=power,
        runs=provenance,
        script_sha256=sha256(__file__),
    )
    args.output.write_text(json.dumps(report, indent=2)+"\n")

    figure, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    masses = [row["mass_min"] for row in rows]
    for index, redshift in enumerate(fine["redshift"]):
        label = f"z = {redshift:g}"
        response = [100*row["missing_response"][index] for row in rows]
        missing = [100*(1-row["resolved_mass"][index]) for row in rows]
        axes[0, 0].semilogx(masses, response, "o-", label=label)
        axes[0, 1].semilogx(masses, missing, "o-", label=label)
        selected = fine["mass"] <= 1e4
        difference = runs["k5"]["sigma"][index]/runs["k9"]["sigma"][index]-1
        axes[1, 0].semilogx(fine["mass"][selected],
                            100*difference[selected], label=label)
        fraction = power["rows"][index]["variance_above_supplied_kmax_fraction"]
        axes[1, 1].semilogx(power["mass"], 100*np.array(fraction),
                            "o-", label=label)
    labels = ["Unresolved bias-weighted response [%]",
              "Ordinary mass budget: 1 − F [%]",
              r"$\sigma$: cutoff $10^5$ versus $10^9$ [%]",
              "Variance from extrapolated power [%]"]
    for axis, label in zip(axes.ravel(), labels):
        axis.set(xlabel=r"Minimum mass [$M_\odot/h$]", ylabel=label)
        axis.grid(alpha=0.2)
    axes[1, 0].set_xlabel(r"Halo mass [$M_\odot/h$]")
    axes[1, 1].set_xlabel(r"Halo mass [$M_\odot/h$]")
    axes[0, 0].legend()
    figure.suptitle("Lower halo-mass cutoff: unchanged fitted halo model\n"
                   "Power beyond the supplied CAMB grid is log-log extrapolated")
    figure.savefig(args.figure, dpi=180)
    figure.savefig(args.figure.with_suffix(".pdf"))
    print(f"Saved microhalo study: {args.output}")


if __name__ == "__main__":
    main()
