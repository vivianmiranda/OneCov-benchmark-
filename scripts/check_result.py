"""Check a small OneCov output and, for shear, its Gaussian normalization.

For an unweighted mean of C_ell in a band containing N integer multipoles,
Gaussian statistics give Var = sum[2*(C_ell+noise)^2/(2*ell+1)]/(fsky*N^2).
This differs from the mode-count weighting in Cocoa's Fourier operator.
The check uses OneCov's saved signal grid and its documented linear
interpolation of the Gaussian numerator, then performs the integer sum
with NumPy. It is an assembly check, not an interpolation-convergence test.

Example: python scripts/check_result.py work/gaussian_shared
"""

import argparse
import configparser
import json
from pathlib import Path

import numpy as np


def main():
    """Save finite/symmetry/positivity and single-source Gaussian checks.

    Returns:
        None. Writes checks.json next to the native result.
    Raises:
        ValueError for incomplete, nonfinite, asymmetric or nonpositive total
        matrices, or a Gaussian normalization mismatch above text precision.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    status = json.loads((run / "run.json").read_text())
    if status["status"] != "completed":
        raise ValueError(f"{run}: run status is {status['status']}, not completed")
    total = np.loadtxt(fname=run / "covariance_matrix.mat")
    if total.ndim != 2 or total.shape[0] != total.shape[1]:
        raise ValueError(f"{run}: expected a square covariance, got {total.shape}")
    if not np.all(np.isfinite(total)) or np.any(np.diag(total) <= 0):
        raise ValueError(f"{run}: covariance must be finite with positive variances")

    # Normalize each entry by the rms of its two measurements. This avoids
    # letting different units or large galaxy variances hide shear errors.
    rms = np.sqrt(np.diag(total))
    correlation = total / rms[:, None] / rms[None, :]
    asymmetry = float(np.max(np.abs(correlation - correlation.T)))
    if asymmetry > 1.e-6:
        raise ValueError(f"{run}: correlation asymmetry {asymmetry:g} exceeds 1e-6")
    minimum_eigenvalue = float(np.linalg.eigvalsh(correlation).min())
    if minimum_eigenvalue <= 0:
        raise ValueError(f"{run}: total is not positive definite; no repair applied")
    checks = {
        "shape": list(total.shape),
        "max_correlation_asymmetry": asymmetry,
        "min_correlation_eigenvalue": minimum_eigenvalue,
        "positive_total": True,
        "components": {},
    }
    component_sum = np.zeros_like(total)
    for filename in sorted(run.glob("covariance_matrix_*.mat")):
        component = np.loadtxt(fname=filename)
        if component.shape != total.shape or not np.all(np.isfinite(component)):
            raise ValueError(f"{filename}: wrong shape or nonfinite covariance")
        checks["components"][filename.name] = {
            "max_absolute_entry": float(np.max(np.abs(component))),
        }
        component_sum += component
    sum_error = (component_sum - total) / rms[:, None] / rms[None, :]
    maximum_sum_error = float(np.max(np.abs(sum_error)))
    checks["component_sum_max_normalized_difference"] = maximum_sum_error
    if maximum_sum_error > 3.e-6:
        raise ValueError("saved G/SSC/cNG components do not sum to the saved total")

    if status["case"] == "shear":
        config = configparser.ConfigParser()
        config.read(filenames=run / "onecov.ini")
        manifest = json.loads((run / "input_manifest.json").read_text())
        grid = config["covELLspace settings"]
        if grid["ell_type_lensing"] != "log":
            raise ValueError("the pilot Gaussian checker currently needs log bands")
        edges = np.geomspace(start=float(grid["ell_min_lensing"]),
                             stop=float(grid["ell_max_lensing"]),
                             num=int(grid["ell_bins_lensing"]) + 1)
        edges = np.unique(edges.astype(int))
        spectrum = np.loadtxt(fname=run / "Cell_kappakappa.ascii")
        ell = spectrum[:, 0]
        signal = spectrum[:, 3]
        fsky = manifest["area_deg2"] * (np.pi / 180.0)**2 / (4.0 * np.pi)

        # Convert the per-arcminute source density into objects per steradian.
        # Noise is sigma_component^2/n_sr, with no extra factor of two.
        arcmin_rad = np.pi / (180.0 * 60.0)
        density_sr = manifest["source_density_arcmin2"] / arcmin_rad**2
        noise = manifest["sigma_e_component"]**2 / density_sr
        numerator = 2.0 * (signal + noise)**2
        expected = np.zeros(shape=total.shape)
        if total.shape != (len(edges)-1, len(edges)-1):
            raise ValueError("matrix dimension does not match one-source band count")

        # Disjoint bands share no Gaussian modes, so only the diagonal is
        # nonzero. Each independent harmonic mode contributes its variance;
        # dividing by the square of the band width forms a mean spectrum.
        for band in range(len(edges)-1):
            modes = np.arange(edges[band], edges[band+1])
            values = np.interp(x=modes, xp=ell, fp=numerator)
            expected[band, band] = np.sum(values / (2*modes+1)) / fsky / len(modes)**2
        gaussian = np.loadtxt(fname=run / "covariance_matrix_gauss.mat")
        scale = np.sqrt(np.diag(expected))
        difference = (gaussian - expected) / scale[:, None] / scale[None, :]
        maximum = float(np.max(np.abs(difference)))
        checks["gaussian_max_normalized_difference"] = maximum
        checks["gaussian_tolerance"] = 3.e-6
        checks["gaussian_scope"] = "uniform-ell bands; saved spectra/text precision"
        if maximum > checks["gaussian_tolerance"]:
            raise ValueError(f"Gaussian assembly discrepancy {maximum:g} exceeds 3e-6")

    (run / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
