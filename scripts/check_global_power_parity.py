"""Check the current power helper against archived native and dense inputs.

No CAMB or covariance calculation is run. The current production helper
refines the saved native power arrays, and every resulting float64 byte is
compared with the dense arrays used by the completed code comparison.
Neither archive nor its original source manifest is modified.
"""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from common import power_table_record, revision, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True,
                        help="original_power.npz from the saved control")
    parser.add_argument("--dense", type=Path, required=True,
                        help="directory with power_tables.npz and manifest.json")
    parser.add_argument("--core", type=Path, default=(
        Path(__file__).resolve().parents[2]
        / "cocoa/Cocoa/external_modules/code/cosmolike_core"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file")

    # Verify the archived files before using them as the expected answers.
    native_report_path = args.native.parent/"report.json"
    native_report = json.loads(native_report_path.read_text())
    manifest_path = args.dense/"manifest.json"
    manifest = json.loads(manifest_path.read_text())
    dense_path = args.dense/"power_tables.npz"
    if sha256(args.native) != native_report["original_power_sha256"]:
        raise ValueError("the archived native power file changed")
    if sha256(dense_path) != manifest["files"]["power_tables.npz"]:
        raise ValueError("the archived dense power file changed")
    with np.load(args.native) as archive:
        native = {name: archive[name] for name in archive.files}
    with np.load(dense_path) as archive:
        dense = {name: archive[name] for name in archive.files}
    refinement = manifest["installed_power"]["power_refinement"]

    # Load just the real preparation helper, without initializing a survey.
    helper_path = args.core/"cosmolike_notebook_utils/covariance/power.py"
    spec = importlib.util.spec_from_file_location("current_power", helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    prepared = helper.refine_power_tables(native, refinement=refinement)
    current = power_table_record(prepared, {"power_refinement": refinement})
    archived = power_table_record(dense, {"power_refinement": refinement})
    if archived != manifest["installed_power"]:
        raise ValueError("dense arrays disagree with their original manifest")

    checks = {}
    for name, expected in archived["arrays"].items():
        checks[name] = current["arrays"][name] == expected
    record = dict(
        status="passed" if all(checks.values()) else "failed",
        scope="Current production helper; archived inputs only; exact bytes",
        arrays_equal=checks, native_nk=len(native["log10k_2D"]),
        installed_power=current, archived_installed_power=archived,
        core=revision(args.core),
        archived_helper_sha256=manifest["power_preparation_sha256"],
        current_helper_sha256=sha256(helper_path),
        native_archive_sha256=sha256(args.native),
        native_report_sha256=sha256(native_report_path),
        dense_archive_sha256=sha256(dense_path),
        dense_manifest_sha256=sha256(manifest_path),
        script_sha256=sha256(__file__),
        limit="Input preparation parity, not a new covariance run or timing")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps({"status": record["status"], "arrays_equal": checks},
                     indent=2))
    if record["status"] != "passed":
        raise SystemExit("power preparation differs; diagnostic preserved")


if __name__ == "__main__":
    main()
