"""Compare shared-spectrum Gaussian band covariances in the Cocoa environment.

This is an assembly test: OneCovariance and Cocoa receive the same LSST
spectra at every integer multipole. A benchmark-only uniform band operator
matches OneCov's estimator; Cocoa's usual mode-count weights are unchanged.
An independent NumPy Wick contraction checks both implementations.

The output contains total, sample-variance, signal-noise and pure-noise
comparisons. It does not validate the supplied spectra, SSC or cNG.
"""

import argparse
import configparser
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

from common import revision, sha256


def load_spectra(inputs, case):
    """Read shared signal into [ell, field, field] with galaxies before shear.

    Arguments:
        inputs = directory containing the exported four-column C_ell files.
        case = 'shear' or the two-lens, one-source '3x2' experiment.
    Returns:
        Integer multipoles, signal cube and observable field pairs.
        The ordering matches OneCov's matrix writer: gg, gm, mm, with bands
        inside each tomographic pair. Cross-lens clustering is retained.
    Raises:
        ValueError if samples are missing, nonfinite or not integer spaced.
    """
    ell = np.loadtxt(fname=inputs / "Cmm.txt")[:, 0]
    if (ell[0] != int(ell[0])
            or not np.array_equal(ell, np.arange(ell[0], ell[-1] + 1))):
        raise ValueError("assembly test requires every integer multipole")
    nlens = 2 if case == "3x2" else 0
    spectra = np.full(shape=(len(ell), nlens + 1, nlens + 1), fill_value=np.nan)
    files = [("Cmm.txt", nlens, nlens)]
    pairs = [(nlens, nlens)]
    if nlens:
        files += [("Cgg.txt", 0, 0), ("Cgm.txt", 0, nlens)]
        pairs = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)]

    # Each file restarts its lens/source IDs at one. Offsets convert them
    # to the common field IDs needed by the crossed Wick contractions.
    for filename, first_offset, second_offset in files:
        table = np.loadtxt(fname=inputs / filename)
        for mode, first, second, value in table:
            node = int(mode - ell[0])
            left = int(first) - 1 + first_offset
            right = int(second) - 1 + second_offset
            spectra[node, left, right] = value
            spectra[node, right, left] = value
    if not np.all(np.isfinite(spectra)):
        raise ValueError("missing or nonfinite crossed spectrum")
    return ell, spectra, np.asarray(pairs, dtype=np.int32)


def wick_reference(spectra, noise, pairs, operators, ell, area_sr):
    """Evaluate the two Gaussian contractions and their noise terms in NumPy.

    Arguments:
        spectra = signal [ell, field, field]; noise = diagonal white powers.
        pairs = measured (A,B) field IDs; operators = band weights [band,ell].
        ell = integer multipoles; area_sr = common footprint in steradians.
    Returns:
        Component matrices, ordered by observable then band.
    """
    nband = len(operators)
    size = len(pairs) * nband
    components = {}
    for name in ("sample_variance", "mixed", "noise"):
        components[name] = np.zeros(shape=(size, size))
    noise_matrix = np.diag(noise)
    mode_count = (2 * ell + 1) * area_sr / (4 * np.pi)

    # Cov(C_AB,C_CD) pairs A with C and B with D, then A with D and B
    # with C. Expanding each (signal + noise) product gives CC, CN and NN.
    # The matrix product below sums their independent harmonic variances
    # with one band weight for each of the two measured spectra.
    for left, (a, b) in enumerate(pairs):
        for right, (c, d) in enumerate(pairs):
            ac, bd = spectra[:, a, c], spectra[:, b, d]
            ad, bc = spectra[:, a, d], spectra[:, b, c]
            nac, nbd = noise_matrix[a, c], noise_matrix[b, d]
            nad, nbc = noise_matrix[a, d], noise_matrix[b, c]
            numerators = {
                "sample_variance": ac * bd + ad * bc,
                "mixed": ac * nbd + nac * bd + ad * nbc + nad * bc,
                "noise": np.full(len(ell), nac * nbd + nad * nbc),
            }
            rows = slice(left * nband, (left + 1) * nband)
            cols = slice(right * nband, (right + 1) * nband)
            for name, numerator in numerators.items():
                block = (operators * (numerator / mode_count)) @ operators.T
                components[name][rows, cols] = block
    components["total"] = sum(components.values())
    return components


def load_onecov_components(run, pairs, nband, nlens):
    """Map OneCov's labeled split-G list to its published matrix ordering.

    Arguments:
        run = completed native output directory; pairs = common field IDs.
        nband = number of bands; nlens = zero for shear, two for 3x2.
    Returns:
        Three component matrices and the higher-precision saved total.
    Raises:
        ValueError for missing entries, unsupported samples or band counts.
    """
    rows = np.loadtxt(fname=run / "covariance_list.dat", dtype=str)
    centers = np.unique(rows[:, 1:3].astype(float))
    if len(centers) != nband:
        raise ValueError("unexpected number of rounded band centers")
    pair_lookup = {tuple(pair): index for index, pair in enumerate(pairs)}
    band_lookup = {center: index for index, center in enumerate(centers)}
    size = len(pairs) * nband
    components = {}
    for name in ("sample_variance", "mixed", "noise"):
        components[name] = np.full(shape=(size, size), fill_value=np.nan)

    # A list row identifies both measured spectra, so no array-layout
    # guess is needed. 'm' denotes this experiment's single shear field.
    for row in rows:
        if row[3] != "1" or row[4] != "1":
            raise ValueError("only one population sample is supported")
        fields = []
        for tracer, bin_id in zip(row[0], row[5:9]):
            field = int(bin_id) - 1
            if tracer == "m":
                field += nlens
            fields.append(field)
        left = pair_lookup[tuple(fields[:2])] * nband
        right = pair_lookup[tuple(fields[2:])] * nband
        left += band_lookup[float(row[1])]
        right += band_lookup[float(row[2])]
        for column, name in enumerate(components, start=10):
            value = float(row[column])
            previous = components[name][left, right]
            # The list can give both orientations of one covariance.
            # Check them before mirroring, rather than hiding asymmetry.
            if np.isfinite(previous) and previous != value:
                raise ValueError("inconsistent repeated Gaussian list entry")
            components[name][left, right] = value
            components[name][right, left] = value
    for matrix in components.values():
        if not np.all(np.isfinite(matrix)):
            raise ValueError("incomplete Gaussian component list")
    components["total"] = np.loadtxt(fname=run / "covariance_matrix_gauss.mat")
    return components


def differences(candidate, reference, rms):
    """Return variance-scaled residuals and fractional nonzero-entry errors.

    Arguments:
        candidate, reference = same-shaped component matrices.
        rms = square root of total reference variances, shared by components.
    Returns:
        Maximum absolute residual in correlation units and maximum fractional
        error where the reference component is nonzero. Zero entries are
        still tested by the first statistic.
    """
    residual = candidate - reference
    nonzero = reference != 0
    fractional = 0.0
    if np.any(nonzero):
        fractional = float(np.max(np.abs(residual[nonzero] / reference[nonzero])))
    return {
        "max_variance_scaled_residual": float(
            np.max(np.abs(residual / rms[:, None] / rms[None, :]))),
        "max_fractional_nonzero_error": fractional,
    }


def main():
    """Compare one native run, preserve matrices and fail on a contract error."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--cocoa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve()
    output = args.output.resolve()
    if output.exists():
        parser.error("choose a new output directory to preserve comparisons")
    status = json.loads((run / "run.json").read_text())
    if (status["status"] != "completed" or status["terms"] != "gaussian"
            or status["spectra"] != "shared-cells"):
        parser.error("requires a completed Gaussian shared-cells run")
    config = configparser.ConfigParser()
    config.read(filenames=run / "onecov.ini")
    inputs = Path(config["tabulated inputs files"]["cell_directory"])
    manifest = json.loads((run / "input_manifest.json").read_text())
    for filename, expected in manifest["files"].items():
        if sha256(filename=inputs / filename) != expected:
            raise ValueError(f"shared input changed: {filename}")

    ell, spectra, pairs = load_spectra(inputs=inputs, case=status["case"])
    grid = config["covELLspace settings"]
    if grid["ell_type_lensing"] != "log":
        raise ValueError("this experiment uses logarithmic output bands")
    for suffix in ("min", "max", "bins", "type"):
        if grid[f"ell_{suffix}_lensing"] != grid[f"ell_{suffix}_clustering"]:
            raise ValueError("this experiment requires identical probe bands")
    edges = np.geomspace(float(grid["ell_min_lensing"]),
                         float(grid["ell_max_lensing"]),
                         int(grid["ell_bins_lensing"]) + 1).astype(int)
    if (np.any(np.diff(edges) <= 0)
            or edges[0] < ell[0]
            or edges[-1] > ell[-1] + 1):
        raise ValueError("empty band or missing integer multipoles")
    operators = np.zeros(shape=(len(edges)-1, len(ell)))
    for band, (lower, upper) in enumerate(zip(edges[:-1], edges[1:])):
        operators[band, (ell >= lower) & (ell < upper)] = 1.0 / (upper - lower)

    arcmin_rad = np.pi / (180 * 60)
    noise = [manifest["sigma_e_component"]**2 * arcmin_rad**2 /
             manifest["source_density_arcmin2"]]
    if status["case"] == "3x2":
        lens_noise = arcmin_rad**2 / np.asarray(manifest["lens_density_arcmin2"])
        noise = list(lens_noise) + noise
    noise = np.asarray(noise)
    area_sr = manifest["area_deg2"] * (np.pi / 180)**2
    reference = wick_reference(spectra=spectra, noise=noise, pairs=pairs,
                               operators=operators, ell=ell, area_sr=area_sr)

    # These production bindings use the same C kernels as survey runs.
    # They accept supplied arrays and need no CAMB or likelihood setup.
    cocoa = args.cocoa.resolve()
    sys.path.insert(0, str(cocoa / "projects" / "lsst_y1"))
    import cosmolike_lsst_y1_interface as ci
    compute = ci.covariance.covariance_gaussian_fourier
    common = dict(pairs=pairs, operators=operators, ell_min=int(ell[0]),
                  area_sr=area_sr)
    result = {
        "total": compute(spectra=spectra, noise=noise, **common),
        "sample_variance": compute(spectra=spectra, noise=np.zeros_like(noise),
                                   **common),
        "noise": compute(spectra=np.zeros_like(spectra), noise=noise, **common),
    }
    result["mixed"] = result["total"] - result["sample_variance"] - result["noise"]
    onecov = load_onecov_components(run=run, pairs=pairs, nband=len(operators),
                                   nlens=len(noise)-1)
    rms = np.sqrt(np.diag(reference["total"]))
    report = {
        "scope": "shared integer-ell spectra; uniform band estimator; Gaussian only",
        "case": status["case"], "shape": list(reference["total"].shape),
        "band_edges": edges.tolist(), "pairs": pairs.tolist(),
        "onecov_run": status, "inputs": manifest,
        "comparer_sha256": sha256(filename=__file__),
        "interface_sha256": sha256(filename=ci.__file__),
        "core": revision(directory=cocoa / "external_modules/code/cosmolike_core"),
        "components": {}, "passed": True,
    }
    for name in reference:
        cocoa_error = differences(candidate=result[name],
                                  reference=reference[name], rms=rms)
        onecov_error = differences(candidate=onecov[name],
                                   reference=reference[name], rms=rms)
        cross_error = differences(candidate=onecov[name],
                                  reference=result[name], rms=rms)
        # OneCov's matrix has seven significant digits; its split list
        # has five. These limits test assembly against output precision.
        tolerance = 1.e-6 if name == "total" else 6.e-5
        precision = ".6e" if name == "total" else ".4e"
        rounded = np.empty_like(reference[name])
        for index in np.ndindex(rounded.shape):
            rounded[index] = float(format(reference[name][index], precision))
        matches_printed = bool(np.array_equal(onecov[name], rounded))
        passed = (cocoa_error["max_variance_scaled_residual"] < 1.e-12
                  and onecov_error["max_variance_scaled_residual"] < tolerance
                  and onecov_error["max_fractional_nonzero_error"] < tolerance)
        report["components"][name] = {
            "cocoa_vs_numpy": cocoa_error, "onecov_vs_numpy": onecov_error,
            "onecov_vs_cocoa": cross_error, "onecov_text_tolerance": tolerance,
            "matches_numpy_at_native_text_precision": matches_printed,
            "passed": passed,
        }
        report["passed"] &= passed

    for label, matrices in (("cocoa", result), ("onecov", onecov)):
        correlation = matrices["total"] / rms[:, None] / rms[None, :]
        minimum = float(np.linalg.eigvalsh(correlation).min())
        report[f"{label}_min_scaled_eigenvalue"] = minimum
        report["passed"] &= minimum > 0
    ratios = eigvalsh(a=onecov["total"], b=result["total"])
    report["generalized_variance_ratio_range"] = [
        float(ratios[0]), float(ratios[-1]),
    ]
    report["max_generalized_variance_change"] = float(
        np.max(np.abs(ratios - 1)))
    report["passed"] = bool(report["passed"])
    output.mkdir(parents=True)
    arrays = {}
    for label, matrices in (("cocoa", result), ("onecov", onecov),
                            ("numpy", reference)):
        for name, matrix in matrices.items():
            arrays[f"{label}_{name}"] = matrix
    np.savez_compressed(file=output / "matrices.npz", **arrays)
    (output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{status['case']}: passed={report['passed']}; largest variance-mode "
          f"change={report['max_generalized_variance_change']:.3g}")
    if not report["passed"]:
        raise SystemExit(f"Comparison failed; inspect {output / 'comparison.json'}")


if __name__ == "__main__":
    main()
