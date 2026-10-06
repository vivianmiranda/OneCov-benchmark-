"""Collect sequential timings and checks; do not run either numerical code.

Projection cases use common resident inputs. Halo moments instead retain
the different native halo prescriptions. Keep these scopes separate when
interpreting a ratio of times.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=Path("work"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file")

    cases = {}
    for case in ("ssc300", "ssc2401", "cng300", "cng601", "halo"):
        row = {}
        matrices = {}
        for backend in ("cocoa", "onecov"):
            records = []
            for suffix in ("", "_repeat"):
                folder = args.work/f"timing_{case}_{backend}{suffix}"
                record = json.loads((folder/"timing.json").read_text())
                if sha256(folder/"outputs.npz") != record["output_sha256"]:
                    raise ValueError(f"timing output changed: {folder}")
                records.append(record)
            samples = [value for record in records
                       for value in record["samples_seconds"]]
            row[backend] = dict(
                mean_seconds=float(np.mean(samples)),
                std_seconds=float(np.std(samples, ddof=1)),
                run_means_seconds=[r["mean_seconds"] for r in records],
                runs=records)
            if case != "halo":
                matrices[backend] = np.load(
                    args.work/f"timing_{case}_{backend}/outputs.npz")["matrix"]
        row["onecov_over_cocoa"] = (
            row["onecov"]["mean_seconds"]/row["cocoa"]["mean_seconds"])
        if case != "halo":
            reference = matrices["onecov"]
            scale = np.sqrt(np.outer(reference.diagonal(), reference.diagonal()))
            difference = np.max(np.abs(matrices["cocoa"]-reference)/scale)
            if difference > 1e-11:
                raise ValueError(f"timed matrices disagree: {case}")
            row["max_variance_scaled_difference"] = float(difference)
        cases[case] = row
        print(case, "; ".join(
            f"{b}: {row[b]['mean_seconds']*1000:.4f} "
            f"+/- {row[b]['std_seconds']*1000:.4f} ms"
            for b in ("cocoa", "onecov")),
            f"; ratio {row['onecov_over_cocoa']:.3f}")

    # The moment benchmark requests the native 200-point power grid.
    # Report convergence both over that whole grid and the physical range
    # used in the preceding halo comparison, rather than hiding its tails.
    wave = np.load(args.work/"halo_onecov_800/inputs.npz")["k"][0]
    left, right = np.triu_indices(len(wave))
    selected = (wave >= 0.001) & (wave <= 10)
    refinements = {}
    for backend, coarse in (("onecov", "timing_halo400_onecov"),
                             ("cocoa", "timing_halo96_cocoa")):
        low = np.load(args.work/coarse/"outputs.npz")
        high = np.load(args.work/f"timing_halo_{backend}/outputs.npz")
        row = {}
        for name in low.files:
            error = np.abs(low[name]/high[name]-1)
            valid = selected if name == "i11" else selected[left]&selected[right]
            row[name] = dict(all_k_max_fractional=float(error.max()),
                            comparison_range_max_fractional=float(
                                error[:, valid].max()))
        refinements[backend] = row

    result = dict(
        scopes=dict(
            ssc="shared-response radial projection, 100x100 shear",
            cng="shared-trispectrum radial projection, 8x8 shear",
            halo="native moments at three redshifts, 200 k nodes; different fits"),
        cases=cases, halo_refinements=refinements,
        halo_k_limits_h_mpc=[float(wave[0]), float(wave[-1])],
        halo_comparison_k_limits_h_mpc=[0.001, 10],
        collector_sha256=sha256(__file__))
    args.output.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
