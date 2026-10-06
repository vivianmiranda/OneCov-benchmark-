"""Isolate numerical FFTLog weighting at unchanged mass and power domains."""

import argparse
import json
from pathlib import Path

import numpy as np

from common import sha256
from report_mass_cutoff import read_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("old", "new", "filter", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output file")
    old, old_record = read_run(args.old)
    new, new_record = read_run(args.new)
    for name in ("boost", "tail_nodes", "tail_panels", "power_table_sha256"):
        if old_record[name] != new_record[name]:
            raise ValueError(f"control differs in {name}")
    for name in ("table_mass_min", "table_mass_max"):
        if old_record["diagnostic_build"][name] != new_record[
                "diagnostic_build"][name]:
            raise ValueError(f"control differs in {name}")
    for name in ("tail_k", "tail_power", "redshift"):
        np.testing.assert_array_equal(old[name], new[name])
    reference = json.loads(args.filter.read_text())
    if reference["source_ingredients_sha256"] != new_record["ingredients_sha256"]:
        raise ValueError("filter does not use the control's exported power")

    # Compare the same physical masses with the native OneCov filter.
    # The old run may contain more exported mass samples than the control.
    errors = {}
    for label, data in (("bias1.5", old), ("bias0.5", new)):
        indices = [np.argmin(abs(data["mass"]/mass-1))
                   for mass in reference["mass"]]
        np.testing.assert_allclose(data["mass"][indices], reference["mass"],
                                   rtol=1e-12, atol=0)
        expected = np.array([row["filter_sigma"][str(float(data["tail_k"][-1]))]
                             for row in reference["rows"]])
        errors[label] = (data["sigma"][:, indices]/expected-1).tolist()
    result = dict(mass=reference["mass"], redshift=new["redshift"].tolist(),
                  fractional_sigma_errors=errors, old_run=old_record,
                  new_run=new_record, filter_sha256=sha256(args.filter),
                  script_sha256=sha256(__file__))
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(errors, indent=2))


if __name__ == "__main__":
    main()
