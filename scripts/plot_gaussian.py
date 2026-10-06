"""Plot the validated Gaussian matrices and variance components.

Run in either prepared environment. Inputs are the reviewed comparison
folders; figures are ordinary PNG/PDF files suitable for the README.
Residuals here use OneCov's saved text, so they expose its output rounding.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


LABELS = [r"$g_1g_1$", r"$g_1g_2$", r"$g_2g_2$",
          r"$g_1\gamma$", r"$g_2\gamma$", r"$\gamma\gamma$"]


def observable_axis(axis, vertical=False):
    """Mark six measured spectra, each containing five adjacent bands."""
    centers = np.arange(6) * 5 + 2
    axis.set_xticks(centers, LABELS)
    for boundary in np.arange(5, 30, 5) - 0.5:
        axis.axvline(boundary, color="0.65", linewidth=0.6)
    if vertical:
        axis.set_yticks(centers, LABELS)
        for boundary in np.arange(5, 30, 5) - 0.5:
            axis.axhline(boundary, color="0.65", linewidth=0.6)


def save_figure(figure, output, name):
    """Save matching screen and vector versions, then release the figure."""
    figure.savefig(output / f"{name}.png", dpi=180, facecolor="white")
    figure.savefig(output / f"{name}.pdf", facecolor="white")
    plt.close(figure)


def main():
    """Create matrix and diagonal-component panels for both multipole ranges."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparisons", type=Path, required=True,
                        help="parent folder containing the reviewed cases")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    cases = [("3x2_30_150", r"$30\leq\ell<150$"),
             ("3x2_1500_1620", r"$1500\leq\ell<1620$")]
    figure, axes = plt.subplots(2, 3, figsize=(12, 8), layout="constrained")
    fractions, fraction_axes = plt.subplots(
        2, 2, figsize=(12, 7), layout="constrained",
        gridspec_kw={"height_ratios": [2, 1]})

    for row, (case, title) in enumerate(cases):
        data = np.load(file=args.comparisons / case / "matrices.npz")
        cocoa = data["cocoa_total"]
        onecov = data["onecov_total"]
        rms = np.sqrt(np.diag(cocoa))
        residual = (onecov - cocoa) / rms[:, None] / rms[None, :]

        # Normalize each code by its own diagonal to show correlations.
        # The difference panel instead uses one common Cocoa normalization,
        # so diagonal variance changes remain visible as well as cross terms.
        for column, (name, matrix) in enumerate(
                (("CoCoA", cocoa), ("OneCovariance", onecov))):
            scale = np.sqrt(np.diag(matrix))
            correlation = matrix / scale[:, None] / scale[None, :]
            image = axes[row, column].imshow(
                correlation, origin="lower", vmin=-1, vmax=1, cmap="RdBu_r")
            axes[row, column].set_title(f"{name}: {title}")
            observable_axis(axis=axes[row, column], vertical=True)
        figure.colorbar(image, ax=list(axes[row, :2]), shrink=0.78,
                        label="Correlation")
        image = axes[row, 2].imshow(residual * 1.e6, origin="lower",
                                    vmin=-0.5, vmax=0.5, cmap="RdBu_r")
        axes[row, 2].set_title("OneCov − CoCoA: saved-text residual")
        observable_axis(axis=axes[row, 2], vertical=True)
        figure.colorbar(image, ax=axes[row, 2], shrink=0.78,
                        label=r"$10^6\Delta C_{ij}/\sqrt{C_{ii}C_{jj}}$")

        # Fractions sum to one along the diagonal. Solid lines are Cocoa;
        # the native OneCov values are markers and should fall on the lines.
        for name, label, color in (
            ("sample_variance", "Sample variance", "#2563a6"),
            ("mixed", "Signal × noise", "#c05621"),
            ("noise", "Pure noise", "#23836c"),
        ):
            values = np.diag(data[f"cocoa_{name}"]) / np.diag(cocoa)
            other = np.diag(data[f"onecov_{name}"]) / np.diag(onecov)
            for start in range(0, 30, 5):
                # Multipoles increase within each spectrum. Do not draw
                # a connecting line across different measured quantities.
                fraction_axes[0, row].plot(
                    np.arange(start, start + 5), values[start:start + 5],
                    color=color, label=label if start == 0 else None)
            fraction_axes[0, row].plot(other, linestyle="none", marker="o",
                                       markersize=3.5, fillstyle="none", color=color)
        fraction_axes[0, row].set_title(title)
        fraction_axes[0, row].set_ylim(-0.03, 1.05)
        fraction_axes[0, row].set_ylabel("Fraction of Gaussian variance")
        fraction_axes[0, row].legend(fontsize=9)
        observable_axis(axis=fraction_axes[0, row])
        fraction_axes[1, row].plot(np.diag(residual) * 1.e6, "o", color="#54278f")
        fraction_axes[1, row].axhline(0, color="0.4", linewidth=0.8)
        fraction_axes[1, row].set_ylabel("Variance difference [ppm]")
        fraction_axes[1, row].set_ylim(-0.5, 0.5)
        observable_axis(axis=fraction_axes[1, row])

    figure.suptitle("LSST Y1 subset · shared spectra · Gaussian covariance")
    fractions.suptitle("Gaussian components · CoCoA lines / OneCov markers")
    save_figure(figure=figure, output=args.output, name="gaussian_matrices")
    save_figure(figure=fractions, output=args.output, name="gaussian_components")

    print(f"Saved Gaussian plots to {args.output}")


if __name__ == "__main__":
    main()
