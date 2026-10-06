"""Separate finite halo-mass limits from multiplicity normalization.

Integrate each code's own fitted functions over peak height. This is a
diagnostic of the extrapolated fit, not a simulation calibration at tiny
halo masses. No bias, mass function or covariance implementation is changed.
"""

import argparse
import contextlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from capture_native import capture_return
from common import revision, sha256


def integrals(nu, multiplicity, bias):
    """Integrate in ln(nu); the measure contributes the extra factor nu."""
    return dict(mass=float(simpson(multiplicity*nu, x=np.log(nu))),
                biased_mass=float(simpson(multiplicity*bias*nu, x=np.log(nu))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("onecov", "cocoa"))
    parser.add_argument("--config", type=Path,
                        default=Path("work/shear_ssc_native/onecov.ini"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nodes", type=int, default=4097)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output directory")
    args.output.mkdir(parents=True)
    root = Path(__file__).resolve().parents[2]
    nu = np.exp(np.linspace(-90, 3.5, args.nodes))
    rows = []
    with (args.output/"native.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            if args.mode == "onecov":
                sys.path.insert(0, str(root/"OneCovariance"))
                from onecov.cov_input import FileInput, Input
                from onecov.cov_polyspectra import PolySpectra
                terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
                    Input().read_input(config_name=str(args.config)))
                tables = FileInput(bias).read_input(config_name=str(args.config))
                model = PolySpectra(0.1, terms, obs, cosmo, bias, hod,
                                    survey, prec, tables)
            else:
                cocoa = root/"cocoa/Cocoa"
                sys.path[:0] = [str(cocoa/"external_modules/code/cosmolike_core"),
                               str(cocoa/"projects/lsst_y1/covariance")]
                import cosmolike_lsst_y1_interface as ci
                from lsst_y1_covariance import configuration, initialize
                initialize(ci, configuration())

            for z in (0.1, 0.5, 1.0):
                if args.mode == "onecov":
                    model.update_mass_func(z, bias, hod, prec)
                    _, captured = capture_return(model.bias, ("nu_new",),
                                                  lambda: model.bias(bias, prec["hm"]))
                    fit = model.mass_func.hmf
                    # hmf supplies the native Tinker multiplicity directly.
                    # Evaluate it at an extended nu grid without constructing
                    # an artificial power spectrum or modifying halo masses.
                    extended = type(fit)(
                        nu2=nu**2, z=z, mass_definition=fit.mass_definition,
                        cosmo=fit.cosmo, delta_c=fit.delta_c, **fit.params)
                    f = extended.fsigma/nu
                    b = model._HaloModel__bias_tinker10_fittfunc(nu)
                    row = dict(z=z, **integrals(nu, f, b),
                               finite_normalization=float(model.norm_bias),
                               finite_nu_range=[float(captured["nu_new"][0]),
                                                float(captured["nu_new"][-1])],
                               bias_multiplier=1/float(model.norm_bias),
                               multiplicity_amplitude=float(extended.normalise))
                else:
                    a = 1/(1+z)
                    f = np.array([ci.fnu(float(n), a) for n in nu])
                    b = np.array([ci.hb1nu(float(n), a) for n in nu])
                    row = dict(z=z, **integrals(nu, f, b))
                rows.append(row)
    repository = (root/"OneCovariance" if args.mode == "onecov" else
                  root/"cocoa/Cocoa/external_modules/code/cosmolike_core")
    record = dict(mode=args.mode, rows=rows, log_nu_limits=[-90, 3.5],
                  nodes=args.nodes, revision=revision(repository),
                  script_sha256=sha256(__file__),
                  config_sha256=sha256(args.config) if args.mode == "onecov" else None)
    (args.output/"result.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
