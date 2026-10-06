"""Collect the completed SSC-followup halo/cNG diagnostic runs.

This records measured comparisons and one-control refinements; it does
not run either covariance engine or assert full-survey convergence.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

from common import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=Path("work"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    def read(directory, filename):
        return json.loads((args.work/directory/filename).read_text())

    def save(filename, record):
        record["collector_sha256"] = sha256(__file__)
        (args.output/filename).write_text(json.dumps(record, indent=2)+"\n")

    bias = {name: read("bias_"+name, "result.json") for name in (
        "onecov", "onecov_fine", "cocoa", "cocoa_fine")}
    save("bias_normalization_20261006.json", bias)

    trispectrum = dict(
        shared_and_native=read("trispectrum_final_256", "comparison.json"),
        coarse_cocoa=read("trispectrum_final_96", "comparison.json"),
        small_corner=read("trispectrum_guard_cocoa", "comparison.json"),
        refinements=[])
    for coarse, fine in (("trispectrum_400", "trispectrum_800"),
                          ("trispectrum_800", "trispectrum_guard_fine"),
                          ("trispectrum_guard_fine", "trispectrum_guard_finer")):
        left = np.load(args.work/coarse/"inputs.npz")["onecov"]
        right = np.load(args.work/fine/"inputs.npz")["onecov"]
        trispectrum["refinements"].append(dict(
            coarse=read(coarse, "export.json"), fine=read(fine, "export.json"),
            max_fractional_by_halo_order=np.max(
                np.abs(left/right-1), axis=(0, 2, 3)).tolist()))
    left = np.load(args.work/"trispectrum_final_96/comparison.npz")["native"]
    right = np.load(args.work/"trispectrum_final_256/comparison.npz")["native"]
    trispectrum["cocoa_96_to_256_max_fractional_by_halo_order"] = np.max(
        np.abs(left/right-1), axis=(0, 2)).tolist()
    save("trispectrum_20261006.json", trispectrum)

    names = ("k9", "k17", "k33", "k65", "k129", "radial", "redshift",
             "redshift_fine", "mass", "corner")
    connected = dict(
        scope="shared-input 8x8 shear cNG; native OneCov 5-band refinements",
        comparisons={name: read("connected_cocoa_"+name, "comparison.json")
                     for name in names}, refinements=[])
    changes = (("k9", "k17"), ("k17", "k33"), ("k33", "k65"),
               ("k65", "k129"), ("k65", "radial"), ("radial", "redshift"),
               ("redshift", "redshift_fine"), ("redshift", "mass"),
               ("redshift", "corner"))
    for coarse, fine in changes:
        left = np.load(args.work/("connected_"+coarse)/"inputs.npz")
        right = np.load(args.work/("connected_"+fine)/"inputs.npz")
        cng = right["band_connected"]
        scale = np.sqrt(np.outer(cng.diagonal(), cng.diagonal()))
        # Generalized eigenvalues inspect all linear combinations of these
        # five bandpowers, rather than only the diagonal variances.
        ratios = eigvalsh(left["band_total"], right["band_total"])
        connected["refinements"].append(dict(
            coarse=coarse, fine=fine,
            max_cng_variance_scaled_change=float(np.max(np.abs(
                left["band_connected"]-cng)/scale)),
            total_variance_ratios=ratios.tolist(),
            max_total_variance_mode_change=float(np.max(np.abs(ratios-1)))))
    save("connected_20261006.json", connected)


if __name__ == "__main__":
    main()
