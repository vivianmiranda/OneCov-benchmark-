"""Summarize repeated full-covariance timings and verify repeat outputs.

Each cutoff uses three fresh processes, including first-use halo tables.
The timing summary excludes initialization and saving; both remain recorded
in the individual run reports. Run this after all numerical jobs finish.
"""

import argparse
import json
from pathlib import Path
import statistics

import numpy as np

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=Path("work"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new summary path")

    summary = dict(
        repeats_per_cutoff=3, threads=8, first_use_tables_included=True,
        initialization_included=False, saving_included=False,
        script_sha256=sha256(__file__), cutoffs={},
    )
    common_power = None

    for exponent in (6, 4, 2):
        prefix = f"cutoff_full_wide{exponent}"
        folders = [args.work/prefix]
        folders += [args.work/f"{prefix}_repeat{n}" for n in (2, 3)]
        records = []

        # A fresh process must reproduce every physical output, not just
        # the diagonal or a few covariance entries. Timing metadata is
        # the only archive field that is expected to differ between runs.
        with np.load(folders[0]/"covariance.npz") as source:
            baseline = {key: source[key] for key in source.files
                        if key != "stages_json"}
        for folder in folders:
            record = json.loads((folder/"report.json").read_text())
            if record["threads"] != 8:
                raise ValueError(f"unexpected thread count: {folder}")
            if sha256(folder/"covariance.npz") != record["covariance_sha256"]:
                raise ValueError(f"covariance hash changed: {folder}")
            if sha256(folder/"power_tables.npz") != record["power_table_sha256"]:
                raise ValueError(f"power hash changed: {folder}")
            if common_power is None:
                common_power = record["power_table_sha256"]
            if record["power_table_sha256"] != common_power:
                raise ValueError(f"power inputs differ: {folder}")
            with np.load(folder/"covariance.npz") as current:
                for key, expected in baseline.items():
                    np.testing.assert_array_equal(current[key], expected)
            records.append(dict(folder=str(folder), **record))

        stages = {}
        for name in records[0]["stages_s"]:
            values = [record["stages_s"][name] for record in records]
            stages[name] = dict(mean_s=statistics.mean(values),
                                sample_stddev_s=statistics.stdev(values),
                                individual_s=values)
        summary["cutoffs"][str(10**exponent)] = dict(
            mass_panels=len(json.loads(str(baseline["settings_json"]))[
                "lnm_edges"])-1,
            repeat_outputs_bitwise_equal=True, stages=stages, runs=records,
        )

    original = summary["cutoffs"]["1000000"]["stages"]["total"]["mean_s"]
    for cutoff in summary["cutoffs"].values():
        cutoff["total_time_ratio_to_1e6"] = (
            cutoff["stages"]["total"]["mean_s"]/original)
    summary["common_power_sha256"] = common_power
    args.output.write_text(json.dumps(summary, indent=2)+"\n")
    for minimum, cutoff in summary["cutoffs"].items():
        print(minimum, cutoff["stages"]["total"],
              "ratio", cutoff["total_time_ratio_to_1e6"])


if __name__ == "__main__":
    main()
