"""Check low-mass sigma coverage using OneCov's native hmf TopHat filter.

The input is CoCoA's saved power including its unchanged edge continuation.
This compares existing integration engines; it introduces no new variance
integrator and does not replace either code's production power spectrum.
Run after sourcing start_onecov.sh.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from hmf.density_field.filters import TopHat

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file")
    report = json.loads((args.input/"report.json").read_text())
    if sha256(args.input/"ingredients.npz") != report["ingredients_sha256"]:
        raise ValueError("input ingredients changed after export")
    data = np.load(args.input/"ingredients.npz")
    mass = np.array([1e2, 1e4, 1e6])
    radius = (3*mass/(4*np.pi*data["rho"]))**(1/3)
    wave = data["tail_k"]
    rows = []
    for row, z in enumerate(data["redshift"]):
        power = data["tail_power"][row]
        values = {}
        for limit in (report["supplied_k_hmpc"][1], 1e5, 1e6, 1e7):
            selected = wave <= limit*(1+1e-14)
            values[str(limit)] = TopHat(wave[selected], power[selected]).sigma(
                radius).tolist()
        coarse = TopHat(wave[::2], power[::2]).sigma(radius)
        fine = np.asarray(values[str(1e7)])
        installed = []
        for m in mass:
            index = np.argmin(np.abs(np.log(data["mass"]/m)))
            if abs(data["mass"][index]/m-1) > 1e-12:
                raise ValueError("export has no matching cutoff mass")
            installed.append(float(data["sigma"][row, index]))
        truncated = np.asarray(values[str(report["supplied_k_hmpc"][1])])
        rows.append(dict(
            redshift=float(z), filter_sigma=values, cocoa_sigma=installed,
            high_k_tail_sigma_fractional=(
                np.asarray(values[str(1e5)])/fine-1).tolist(),
            filter_refinement_fractional=(coarse/fine-1).tolist(),
            cocoa_vs_shared_filter_fractional=(
                np.asarray(installed)/fine-1).tolist(),
            variance_above_supplied_kmax_fraction=(
                1-(truncated/fine)**2).tolist(),
        ))
    result = dict(
        mass=mass.tolist(), radius_mpc_h=radius.tolist(), rows=rows,
        supplied_kmax_hmpc=report["supplied_k_hmpc"][1],
        source_report_sha256=sha256(args.input/"report.json"),
        source_ingredients_sha256=sha256(args.input/"ingredients.npz"),
        script_sha256=sha256(__file__),
    )
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
