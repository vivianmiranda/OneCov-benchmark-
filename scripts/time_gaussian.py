"""Time Gaussian component assembly from the validated resident spectra.

Run separately in the Cocoa and OneCov environments, never concurrently.
Both routes calculate sample variance, mixed signal/noise and pure noise.
Imports, input reading, cosmology/halo setup and output writing are outside
the assembly timer. The first call is recorded separately from repeats.
This is a small shared-spectrum benchmark, not a full-survey CLI timing.
"""

import argparse
import configparser
import contextlib
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np

from common import revision, sha256
from compare_gaussian import differences, load_spectra


def cocoa_call(inputs, report):
    """Prepare the production binding and return a callable for all G parts.

    Inputs are the validated export directory and comparison report.
    The callable includes allocations and three C kernel calls. Subtracting
    CC and NN from the full (C+N) contraction isolates the mixed component.
    """
    import cosmolike_lsst_y1_interface as ci

    ell, spectra, pairs = load_spectra(inputs=inputs, case=report["case"])
    manifest = report["inputs"]
    arcmin_rad = np.pi / (180 * 60)
    noise = [manifest["sigma_e_component"]**2 * arcmin_rad**2 /
             manifest["source_density_arcmin2"]]
    if report["case"] == "3x2":
        lens_noise = arcmin_rad**2 / np.asarray(manifest["lens_density_arcmin2"])
        noise = list(lens_noise) + noise
    noise = np.asarray(noise)
    edges = report["band_edges"]
    weights = np.zeros(shape=(len(edges)-1, len(ell)))
    for band, (lower, upper) in enumerate(zip(edges[:-1], edges[1:])):
        weights[band, (ell >= lower) & (ell < upper)] = 1.0 / (upper - lower)
    zero_noise = np.zeros_like(noise)
    zero_spectra = np.zeros_like(spectra)
    arguments = dict(pairs=pairs, operators=weights, ell_min=int(ell[0]),
                     area_sr=manifest["area_deg2"] * (np.pi / 180)**2)
    compute = ci.covariance.covariance_gaussian_fourier

    def assemble():
        """Recompute three contractions; reuse inputs, never cached outputs."""
        total = compute(spectra=spectra, noise=noise, **arguments)
        sample = compute(spectra=spectra, noise=zero_noise, **arguments)
        pure_noise = compute(spectra=zero_spectra, noise=noise, **arguments)
        return dict(sample_variance=sample, noise=pure_noise, total=total,
                    mixed=total - sample - pure_noise)

    return assemble, {"interface_sha256": sha256(filename=ci.__file__)}


def onecov_call(config, onecov):
    """Initialize unmodified OneCov outside the timer and use its public G API.

    Returns the assembly callable and initialized object. Its native output
    is 18 blocks: three components for each of six tracer combinations.
    Formatting those blocks as a single matrix is checked outside the timer.
    """
    sys.path.insert(0, str(onecov))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_ell_space import CovELLSpace

    terms, obs, output, cosmo, bias, ia, hod, survey, prec = Input().read_input(
        config_name=str(config))
    tables = FileInput(bias).read_input(config_name=str(config))
    covariance = CovELLSpace(terms, obs, output, cosmo, bias, ia, hod,
                             survey, prec, tables)

    def assemble():
        """Recompute all native Gaussian components and their band integrals."""
        return covariance.covELL_gaussian(
            covELLspacesettings=obs["ELLspace"],
            survey_params_dict=survey, calc_prefac=True)

    return assemble, covariance


def onecov_matrices(blocks, pairs, nband, nlens):
    """Reorder native output for validation only; perform no new integrals.

    Each native block has axes [band,band,sample,sample,A,B,C,D]. There is
    one population sample. The three output matrices use the already
    validated observable order; reversing the two spectra transposes a block.
    """
    groups = {"gggg": 0, "gggm": 1, "ggmm": 2,
              "gmgm": 3, "mmgm": 4, "mmmm": 5}
    components = {}
    size = len(pairs) * nband
    for role, name in enumerate(("sample_variance", "mixed", "noise")):
        matrix = np.zeros(shape=(size, size))
        for left, left_pair in enumerate(pairs):
            for right, right_pair in enumerate(pairs):
                fields = list(left_pair) + list(right_pair)
                kinds = "".join("g" if field < nlens else "m"
                                for field in fields)
                transpose = kinds not in groups
                if transpose:
                    fields = fields[2:] + fields[:2]
                    kinds = kinds[2:] + kinds[:2]
                block = blocks[3 * groups[kinds] + role]
                if np.isscalar(block):
                    continue  # physically absent noise cross-contraction
                bins = tuple(field if field < nlens else field - nlens
                             for field in fields)
                values = block[(slice(None), slice(None), 0, 0) + bins]
                if transpose:
                    values = values.T
                rows = slice(left * nband, (left + 1) * nband)
                cols = slice(right * nband, (right + 1) * nband)
                matrix[rows, cols] = values
        components[name] = matrix
    components["total"] = sum(components.values())
    return components


def main():
    """Record repeated assembly times, then verify the in-memory components."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison", type=Path)
    parser.add_argument("--backend", choices=("cocoa", "onecov"), required=True)
    parser.add_argument("--onecov", type=Path, default=Path("../OneCovariance"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=31)
    parser.add_argument("--batch-size", type=int, default=None)
    args = parser.parse_args()
    if args.output.exists() or args.repeats < 2:
        parser.error("choose a new output directory and at least two repeats")
    batch = args.batch_size
    if batch is None:
        batch = 100 if args.backend == "cocoa" else 10
    if batch < 1 or int(os.environ.get("OMP_NUM_THREADS", "0")) < 1:
        parser.error("positive batch size and OMP_NUM_THREADS are required")
    output = args.output.resolve()
    comparison = args.comparison.resolve()
    report = json.loads((comparison / "comparison.json").read_text())
    if not report["passed"]:
        parser.error("only a validated shared-spectrum comparison can be timed")
    native_config = Path(report["onecov_run"]["command"][-1])
    config = configparser.ConfigParser()
    config.read(filenames=native_config)
    inputs = Path(config["tabulated inputs files"]["cell_directory"])
    for filename, expected in report["inputs"]["files"].items():
        if sha256(filename=inputs / filename) != expected:
            raise ValueError(f"shared input changed: {filename}")
    output.mkdir(parents=True)
    # The parser/constructor may save resolved settings and spectra. Send
    # those files to this new directory, preserving the original run.
    config["output settings"]["directory"] = str(output)
    config_path = output / "onecov.ini"
    with config_path.open(mode="w") as stream:
        config.write(stream)

    with (output / "initialization.log").open(mode="w") as log:
        with contextlib.redirect_stdout(log):
            started = time.perf_counter()
            if args.backend == "cocoa":
                assemble, provenance = cocoa_call(inputs=inputs, report=report)
            else:
                assemble, covariance = onecov_call(
                    config=config_path, onecov=args.onecov.resolve())
                provenance = revision(directory=args.onecov)
            setup_seconds = time.perf_counter() - started
            started = time.perf_counter()
            first = assemble()
            first_seconds = time.perf_counter() - started

    # Batches reduce clock/Python-loop overhead for these small matrices.
    # Every call still allocates and computes its own outputs. Record the
    # first call above, rather than silently hiding it as a warm-up.
    samples = []
    with open(os.devnull, mode="w") as quiet:
        with contextlib.redirect_stdout(quiet):
            for repeat in range(args.repeats):
                started = time.perf_counter()
                for call in range(batch):
                    result = assemble()
                    del result
                samples.append((time.perf_counter() - started) / batch)

    if args.backend == "onecov":
        first = onecov_matrices(blocks=first, pairs=report["pairs"],
                                nband=len(report["band_edges"])-1,
                                nlens=2 if report["case"] == "3x2" else 0)
    reference = np.load(file=comparison / "matrices.npz")
    rms = np.sqrt(np.diag(reference["numpy_total"]))
    checks = {}
    for name, matrix in first.items():
        checks[name] = differences(candidate=matrix,
                                   reference=reference[f"numpy_{name}"], rms=rms)
        if checks[name]["max_variance_scaled_residual"] > 1.e-12:
            raise ValueError(f"{name}: timed calculation differs from reference")
    timing = {
        "backend": args.backend, "case": report["case"],
        "band_edges": report["band_edges"], "shape": report["shape"],
        "scope": "resident shared spectra to three Gaussian components; allocations included; OneCov output repacking excluded",
        "setup_scope": "imports inside backend setup, input parsing and initialization; not a process-start timer",
        "setup_seconds": setup_seconds, "first_call_seconds": first_seconds,
        "samples_seconds": samples, "batch_size": batch,
        "mean_seconds": float(np.mean(samples)),
        "std_seconds": float(np.std(samples, ddof=1)),
        "median_seconds": float(np.median(samples)),
        "omp_threads": int(os.environ["OMP_NUM_THREADS"]),
        "thread_environment": {}, "platform": platform.platform(),
        "implementation": provenance, "checks": checks,
        "comparison_sha256": sha256(filename=comparison / "comparison.json"),
        "script_sha256": sha256(filename=__file__),
    }
    for name in ("OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC",
                 "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        timing["thread_environment"][name] = os.environ.get(name)
    np.savez_compressed(file=output / "components.npz", **first)
    (output / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
    print(f"{args.backend} {report['case']}: "
          f"{timing['mean_seconds']*1000:.4f} ms per assembly; checks passed")


if __name__ == "__main__":
    main()
