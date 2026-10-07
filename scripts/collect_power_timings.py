"""Verify and collect the global-power timing refresh without running a code.

Keep resident-input component times separate from first-use complete pilots.
Every saved matrix and input must pass its checks before a report is written.
"""

import argparse
import configparser
import json
from pathlib import Path

import numpy as np

from common import sha256


CASES = (
    "gaussian_shear_30_150", "gaussian_3x2_30_150",
    "gaussian_3x2_1500_1620", "ssc_601", "ssc_1201",
    "connected_300", "connected_601", "halo_800",
)
COMPLETE = ("fourier_cocoa", "fourier_onecov", "real_cocoa", "real_onecov")


def read_json(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def arrays(path):
    with np.load(path) as saved:
        result = {name: saved[name] for name in saved.files}
    require(all(np.isfinite(x).all() for x in result.values()),
            f"Nonfinite output: {path}")
    return result


def equal_arrays(left, right):
    require(left.keys() == right.keys(), "Array names differ")
    equal = {name: np.array_equal(value, right[name])
             for name, value in left.items()}
    require(all(equal.values()), f"Arrays differ: {equal}")
    return equal


def ini_values(path):
    config = configparser.ConfigParser()
    require(config.read(path), f"Missing configuration: {path}")
    return {section: dict(config[section]) for section in config.sections()}


def same_ini(left, right):
    # Only the output destination changes between a timed run and its
    # accuracy archive. All physical and numerical settings must survive.
    a, b = ini_values(left), ini_values(right)
    a["output settings"].pop("directory")
    b["output settings"].pop("directory")
    require(a == b, f"Native settings changed: {left}")


def max_scaled_difference(left, right, diagonal):
    rms = np.sqrt(diagonal)
    require(np.all(rms > 0), "Nonpositive reference variance")
    return float(np.max(np.abs(left-right)/rms[:, None]/rms[None, :]))


def component_record(stage, root, installed):
    command = stage["command"]
    folder = root/stage["name"]
    record = read_json(folder/"timing.json")
    require(record["script_sha256"] == stage["script_sha256"],
            "Timing script fingerprint changed")
    require(sha256(command[1]) == stage["script_sha256"],
            "Timing script has changed since execution")
    environment = record["thread_environment"]
    require(int(record.get("omp_threads",
                           environment.get("OMP_NUM_THREADS"))) == 8,
            "Unexpected OpenMP count")
    require(environment["OPENBLAS_NUM_THREADS"] == "1",
            "Unexpected BLAS count")

    if stage["name"].startswith("gaussian"):
        source = Path(command[2])
        comparison = read_json(source/"comparison.json")
        require(comparison["passed"], "Unvalidated Gaussian input")
        require(record["comparison_sha256"] == sha256(source/"comparison.json"),
                "Gaussian comparison fingerprint changed")
        native = Path(comparison["onecov_run"]["command"][-1])
        same_ini(folder/"onecov.ini", native)
        inputs = Path(ini_values(native)["tabulated inputs files"]
                      ["cell_directory"])
        for name, expected in comparison["inputs"]["files"].items():
            require(sha256(inputs/name) == expected, f"Input changed: {name}")
        require(all(c["max_variance_scaled_residual"] <= 1e-12
                    for c in record["checks"].values()),
                "Gaussian numerical check failed")
        output = folder/"components.npz"
        fingerprint = record["comparison_sha256"]
    else:
        source = Path(command[3])
        export = read_json(source/"export.json")
        require(record["input_export"] == export, "Export metadata changed")
        require(sha256(source/"inputs.npz") == export["inputs_sha256"],
                "Supplied projection/halo inputs changed")
        require(all(np.isfinite(v) and v <= 1e-11
                    for v in record["checks"].values()),
                "Projection/halo numerical check failed")
        output = folder/"outputs.npz"
        require(sha256(output) == record["output_sha256"],
                "Timing output fingerprint changed")
        fingerprint = export["inputs_sha256"]
        if record["stage"] == "halo" and record["backend"] == "cocoa":
            require(record["implementation"]["installed_power"] == installed,
                    "Halo power differs from the shared export")
        if record["stage"] == "halo" and record["backend"] == "onecov":
            # The ingredient export requested no full trispectrum. Enable
            # its parser and select the native 200-point moment grid, as
            # the timer does; retain every other input and numerical choice.
            expected = ini_values(source/"onecov.ini")
            expected["covariance terms"]["nongauss"] = "True"
            wave = arrays(source/"inputs.npz")["k"][0]
            expected["trispec evaluation"].update(
                log10k_bins=str(len(wave)),
                log10k_min=str(np.log10(wave[0])),
                log10k_max=str(np.log10(wave[-1])),
                small_k_damping_for1h="none")
            actual = ini_values(folder/"onecov.ini")
            for config in (expected, actual):
                config["output settings"].pop("directory")
            require(actual == expected, "Native halo settings changed")
    return record, arrays(output), fingerprint, sha256(output)


def collect_components(stages, root, installed):
    cases = {}
    for case in CASES:
        row, values, fingerprints = {}, {}, []
        for backend in ("cocoa", "onecov"):
            records, outputs = [], []
            for repeat in (1, 2):
                name = f"{case}_{backend}_{repeat}"
                record, output, fingerprint, output_sha = component_record(
                    stages[name], root, installed)
                records.append(dict(
                    record=record, output_sha256=output_sha,
                    timing_sha256=sha256(root/name/"timing.json")))
                outputs.append(output)
                fingerprints.append(fingerprint)
            equal_arrays(*outputs)
            values[backend] = outputs[0]

            # Samples already contain seconds per call, after division by
            # the batch size. Pool both fresh processes with equal weight.
            samples = np.concatenate([r["record"]["samples_seconds"]
                                      for r in records])
            require(len(samples) == (22 if case == "halo_800" else 62),
                    "Wrong number of repeated timing samples")
            require(np.all(np.isfinite(samples) & (samples > 0)),
                    "Invalid timing sample")
            row[backend] = dict(
                mean_seconds=float(samples.mean()),
                std_seconds=float(samples.std(ddof=1)),
                process_means_seconds=[r["record"]["mean_seconds"]
                                       for r in records],
                samples=len(samples), repeats_bitwise_equal=True, runs=records)
        require(len(set(fingerprints)) == 1, "Different inputs between runs")
        row["onecov_over_cocoa"] = (
            row["onecov"]["mean_seconds"]/row["cocoa"]["mean_seconds"])
        if case != "halo_800":
            reference = values["onecov"]
            diagonal = reference["total" if case.startswith("gaussian")
                                 else "matrix"].diagonal()
            differences = {
                name: max_scaled_difference(value, reference[name], diagonal)
                for name, value in values["cocoa"].items()}
            require(max(differences.values()) <= 1e-11,
                    f"Cross-code supplied-input mismatch: {case}")
            row["cross_code_max_variance_scaled_difference"] = differences
        cases[case] = row
    return cases


def collect_complete(root, accuracy, installed):
    complete = {}
    # Revision records and elapsed times may change. All numerical settings,
    # input hashes, output hashes and saved diagnostics must remain equal.
    allowed = {"core", "project", "setup_seconds", "construction_seconds",
               "stages_seconds", "config_sha256", "timing_note"}
    for name in COMPLETE:
        folder, archive = root/name, accuracy/name
        report = read_json(folder/"report.json")
        original = read_json(archive/"report.json")
        require(report.keys() ^ original.keys() <= allowed,
                "Unexpected report schema change")
        different = [key for key in report
                     if report[key] != original.get(key)]
        require(set(different) <= allowed, f"Settings changed: {different}")
        for location, metadata in ((folder, report), (archive, original)):
            require(sha256(location/"covariance.npz") ==
                    metadata["covariance_sha256"], "Covariance hash changed")
        equal = equal_arrays(arrays(folder/"covariance.npz"),
                             arrays(archive/"covariance.npz"))
        require(report["minimum_total_correlation_eigenvalue"] > 0,
                "Total positivity check failed")
        if name.endswith("cocoa"):
            require(report["installed_power"] == installed,
                    "Complete power differs from the shared export")
        else:
            same_ini(folder/"onecov.ini", archive/"onecov.ini")
        complete[name] = dict(
            report=report, accuracy_report=original, arrays_bitwise_equal=equal,
            total_seconds=(report["setup_seconds"]
                           + report["construction_seconds"]),
            report_sha256=sha256(folder/"report.json"),
            accuracy_report_sha256=sha256(archive/"report.json"))
    return complete


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timings", type=Path, required=True)
    parser.add_argument("--accuracy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file")
    ledger = read_json(args.timings/"timing_run.json")
    require(ledger["status"] == "completed", "Timing runner did not complete")
    require(len(ledger["stages"]) == 36, "Expected 36 timing stages")
    require(all(s["status"] == "completed" and s["returncode"] == 0
                for s in ledger["stages"]), "A timing stage failed")
    stages = {s["name"]: s for s in ledger["stages"]}
    manifest = read_json(args.accuracy/"inputs/manifest.json")
    for name, expected in manifest["files"].items():
        require(sha256(args.accuracy/"inputs"/name) == expected,
                f"Shared manifest input changed: {name}")
    installed = manifest["installed_power"]
    require(installed["nk"] == 11993, "Expected the adopted dense power grid")
    cases = collect_components(stages, args.timings, installed)
    complete = collect_complete(args.timings, args.accuracy, installed)
    result = dict(
        status="verified", hardware="Apple M2 Pro, macOS 13.7.5",
        scopes={
            "gaussian": "Resident spectra; three Gaussian components",
            "ssc": "Shared response radial projection, 100x100 shear",
            "connected": "Shared trispectrum radial projection, 8x8 shear",
            "halo": "Native moments, 3 redshifts, 200 k; different halo fits",
            "complete": "Single fresh process per code and space; setup plus "
                        "first-use construction, excluding plotting/writing",
        },
        cases=cases, complete=complete, timing_ledger=ledger,
        input_manifest_sha256=sha256(args.accuracy/"inputs/manifest.json"),
        collector_sha256=sha256(__file__))
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    for name, row in cases.items():
        print(name, [(b, row[b]["mean_seconds"]*1000,
                      row[b]["std_seconds"]*1000) for b in ("cocoa", "onecov")],
              "ratio", row["onecov_over_cocoa"])
    for name, row in complete.items():
        print(name, "setup + construction =", row["total_seconds"], "s")


if __name__ == "__main__":
    main()
