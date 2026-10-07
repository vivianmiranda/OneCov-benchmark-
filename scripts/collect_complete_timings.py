"""Verify a four-pilot timing rerun against accuracy and earlier timing runs.

This reads saved results only. It neither computes covariance matrices nor
replaces the earlier timing record.
"""

import argparse
import json
from pathlib import Path

from collect_power_timings import COMPLETE, collect_complete, read_json, require
from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timings", type=Path, required=True)
    parser.add_argument("--accuracy", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file; preserve earlier records")

    ledger = read_json(args.timings/"timing_run.json")
    require(ledger["status"] == "completed", "Timing runner is unfinished")
    require([s["name"] for s in ledger["stages"]] == list(COMPLETE),
            "Expected exactly the four complete pilot runs")
    for stage in ledger["stages"]:
        require(stage["status"] == "completed" and stage["returncode"] == 0,
                f"Timing stage failed: {stage['name']}")
        require(sha256(stage["command"][1]) == stage["script_sha256"],
                "Timing script changed since execution")

    # Check every component, observable coordinate and signal twice: first
    # against the accuracy archive, then against the earlier timed output.
    # Native settings and input fingerprints are checked by the same helper.
    manifest = read_json(args.accuracy/"inputs/manifest.json")
    installed = manifest["installed_power"]
    complete = collect_complete(args.timings, args.accuracy, installed)
    previous = collect_complete(args.timings, args.previous, installed)
    changes = {}
    for name in COMPLETE:
        new = complete[name]["report"]
        old = previous[name]["accuracy_report"]
        changes[name] = {
            quantity: dict(previous_seconds=old[quantity],
                           current_seconds=new[quantity],
                           percent_change=100*(new[quantity]/old[quantity]-1))
            for quantity in ("setup_seconds", "construction_seconds")}
        previous[name] = dict(
            report=old, report_sha256=previous[name]["accuracy_report_sha256"],
            arrays_bitwise_equal=previous[name]["arrays_bitwise_equal"])

    state = read_json(args.state)
    result = dict(
        status="verified", scope="Four complete source-bin pilot timings; "
        "one fresh process per code and space, no component timing rerun",
        fan_condition=state["fan_condition"],
        quiet_cpu_idle_percent=state["quiet_cpu_idle_percent"],
        other_numerical_compiler_jobs=state["other_numerical_compiler_jobs"],
        complete=complete, previous=previous, timing_changes=changes,
        timing_ledger=ledger, collector_sha256=sha256(__file__),
        verification_helper_sha256=sha256(
            Path(__file__).with_name("collect_power_timings.py")))
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    for name, row in changes.items():
        value = row["construction_seconds"]
        print(name, f"{value['current_seconds']:.4f} s; "
              f"{value['percent_change']:+.2f}% versus previous")


if __name__ == "__main__":
    main()
