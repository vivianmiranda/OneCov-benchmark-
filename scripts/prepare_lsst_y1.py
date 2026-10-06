"""Export an LSST Y1 source and two lens bins for bounded OneCov pilots.

Run inside the activated Cocoa environment. The project supplies its
cosmology, catalog densities and redshift distributions. Only the selected
tomographic columns are exported; their densities are not reassigned.

The shared C_ell files contain signal alone, from Cocoa's production
covariance interface. They isolate Gaussian assembly from differences in
how the codes project matter power. The CAMB power file supports a later,
separate test of that projection. Neither export is a full covariance run.

Example from the benchmark root, after starting Cocoa:
    python scripts/prepare_lsst_y1.py --cocoa ../cocoa/Cocoa \
        --output work/lsst_y1
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from common import revision, sha256


def export_distribution(filename, bin_numbers, output):
    """Copy selected LSST n(z) columns without changing their z convention.

    Arguments:
        filename = LSST text file with midpoint z followed by five n(z) columns.
        bin_numbers = original 1-based tomographic bins, from 1 through 5.
        output = destination file; columns are z and the selected n(z).
    Returns:
        Original redshift samples, used also for the constant bias table.
    Raises:
        ValueError for malformed, negative or unordered distribution data.
    """
    data = np.loadtxt(fname=filename)
    if data.ndim != 2 or data.shape[1] != 6:
        raise ValueError(f"{filename}: expected z plus five LSST n(z) columns")
    if not np.all(np.isfinite(data)) or np.any(data < 0):
        raise ValueError(f"{filename}: z and densities must be finite and nonnegative")
    if (np.any(np.diff(data[:, 0]) <= 0)
            or np.any(np.max(data[:, bin_numbers], axis=0) <= 0)):
        raise ValueError(f"{filename}: need increasing z and a nonempty bin")
    np.savetxt(fname=output, X=data[:, [0] + bin_numbers], header="z_mid n_z")
    return data[:, 0]


def main():
    """Write shared LSST inputs and a manifest into a new output directory.

    Returns:
        None. Writes text spectra/distributions and input provenance.
    Raises:
        Errors from missing Cocoa dependencies or an existing destination.
        Activate Cocoa and choose a new output directory before retrying.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cocoa", type=Path, required=True,
                        help="runtime root containing projects/lsst_y1")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-bin", type=int, choices=range(1, 6), default=3)
    parser.add_argument("--lens-bins", type=int, nargs=2, choices=range(1, 6),
                        default=[1, 2], help="two distinct original LSST lens bins")
    parser.add_argument("--integer-ell-range", type=int, nargs=2,
                        metavar=("FIRST", "LAST"),
                        help="export every integer ell, inclusive; assembly test")
    args = parser.parse_args()
    cocoa = args.cocoa.resolve()
    output = args.output.resolve()
    project = cocoa / "projects" / "lsst_y1"
    core = cocoa / "external_modules" / "code" / "cosmolike_core"
    if not (project / "covariance" / "lsst_y1_covariance.py").is_file():
        parser.error(f"{project}: missing LSST covariance setup; check --cocoa")
    if output.exists():
        parser.error(f"{output} exists; choose a new --output to preserve inputs")
    if args.lens_bins[0] == args.lens_bins[1]:
        parser.error("--lens-bins must select two distinct LSST populations")
    if args.integer_ell_range is not None:
        first, last = args.integer_ell_range
        if first < 2 or last <= first:
            parser.error("--integer-ell-range requires 2 <= FIRST < LAST")

    # Project initialization installs exactly its forecast cosmology in C.
    # Its production bindings avoid the notebook Armadillo array conversion.
    sys.path.insert(0, str(project / "covariance"))
    sys.path.insert(0, str(core))
    import camb
    import cosmolike_lsst_y1_interface as ci
    from lsst_y1_covariance import configuration, initialize
    from cosmolike_notebook_utils.covariance.gaussian import (
        limber_spectra,
        observed_spectra,
    )

    settings = configuration(gaussian={"nonlimber": False, "ia": "none"})
    tables = initialize(interface=ci, settings=settings)
    cosmology = settings["cosmology"]

    # sigma8 is the rms linear density contrast in an 8 Mpc/h sphere at z=0.
    # OneCov uses it instead of As. Obtain it from CAMB with the project's
    # As and background, rather than substituting a rounded amplitude.
    h = cosmology["H0"] / 100.0
    pars = camb.set_params(
        H0=cosmology["H0"],
        ombh2=cosmology["omegab"] * h**2,
        omch2=(cosmology["omegam"] - cosmology["omegab"]) * h**2,
        mnu=0.0,
        omk=0.0,
        tau=0.06,
        As=cosmology["As_1e9"] * 1.e-9,
        ns=cosmology["ns"],
        w=cosmology["w"],
        wa=cosmology["w0pwa"] - cosmology["w"],
        lmax=10,
    )
    pars.set_matter_power(redshifts=[0.0], kmax=20.0)
    camb_result = camb.get_results(params=pars)
    sigma8 = float(camb_result.get_sigma8_0())

    # A dense spectrum grid spans the pilot's five output bands. A source
    # leg needs the shear spin factor; noise is added later by OneCov.
    ell = np.geomspace(start=10.0, stop=5000.0, num=257)
    if args.integer_ell_range is not None:
        # Gaussian band sums run over integer multipoles. Supplying each
        # one removes spectrum interpolation from an assembly comparison.
        first, last = args.integer_ell_range
        ell = np.arange(first, last + 1, dtype=float)
    snapshot = limber_spectra(
        interface=ci.covariance,
        ell=ell,
        a_edges=settings["a_edges"],
        nquad=settings["radial_nquad"],
        nwindow=settings["nwindow"],
    )
    spectra = observed_spectra(spectra=snapshot["spectra"], ell=ell, nlens=5)
    lenses = np.array(args.lens_bins) - 1
    source = 5 + args.source_bin - 1

    output.mkdir(parents=True)
    lens_file = project / settings["lens_file"]
    source_file = project / settings["source_file"]
    z_lens = export_distribution(
        filename=lens_file, bin_numbers=args.lens_bins, output=output / "lens.nz",
    )
    export_distribution(
        filename=source_file, bin_numbers=[args.source_bin],
        output=output / "source.nz",
    )
    bias = np.asarray(settings["bias"])[lenses]
    bias_values = np.tile(bias, (len(z_lens), 1))
    np.savetxt(fname=output / "bias.txt", X=np.column_stack((z_lens, bias_values)),
               header="z bias")

    # Every ordered pair is required, including crossed galaxy bins which
    # might not be in a data vector. Local IDs restart at 1 for each tracer;
    # the manifest retains the corresponding original LSST IDs.
    for name, first_fields, second_fields in (
        ("Cmm.txt", [source], [source]),
        ("Cgm.txt", lenses, [source]),
        ("Cgg.txt", lenses, lenses),
    ):
        rows = []
        for index, mode in enumerate(ell):
            for first_id, first in enumerate(first_fields, start=1):
                for second_id, second in enumerate(second_fields, start=1):
                    rows.append([mode, first_id, second_id,
                                 spectra[index, first, second]])
        values = np.asarray(rows)
        np.savetxt(fname=output / name, X=values, header="ell tomo_i tomo_j C_ell")

    # CAMB's table uses k in h/Mpc and P in (Mpc/h)^3. Its flat storage
    # advances redshift first; reshape before writing OneCov's z-major file.
    k = 10.0**tables["log10k_2D"]
    z = tables["z_2D"]
    power = np.exp(tables["lnP_nonlinear"].reshape((len(z), len(k)), order="F"))
    # Keep the first CAMB node above the survey's maximum redshift too.
    # Its value closes the interpolation interval; discarding it would
    # leave the final part of n(z) outside the supplied power table.
    z_max = 1.0 / settings["a_edges"][0] - 1.0
    last = int(np.searchsorted(z, z_max, side="left")) + 1
    selected = np.arange(len(z)) < last
    power_columns = np.column_stack((np.repeat(z[selected], len(k)),
                                     np.tile(k, np.count_nonzero(selected)),
                                     power[selected].ravel()))
    np.savetxt(fname=output / "Pmm.txt", X=power_columns,
               header="z k_h_Mpc P_Mpc_over_h_cubed")

    manifest = {
        "survey": "LSST Y1 forecast subset",
        "source_bin_1based": args.source_bin,
        "lens_bins_1based": args.lens_bins,
        "area_deg2": settings["area_deg2"],
        "source_density_arcmin2": settings["source_density_arcmin2"][args.source_bin-1],
        "lens_density_arcmin2": np.asarray(
            settings["lens_density_arcmin2"])[lenses].tolist(),
        "sigma_e_component": settings["sigma_e_component"][args.source_bin-1],
        "bias": bias.tolist(),
        "cosmology": cosmology,
        "sigma8": sigma8,
        "camb_version": camb.__version__,
        "camb_source": revision(directory=cocoa / "external_modules" /
                                "code" / "CAMB"),
        "gaussian": settings["gaussian"],
        "ell_grid": {
            "first": float(ell[0]), "last": float(ell[-1]),
            "count": len(ell),
            "spacing": "integer" if args.integer_ell_range else "log",
        },
        "radial_nquad": settings["radial_nquad"],
        "nwindow": settings["nwindow"],
        "a_edges": settings["a_edges"].tolist(),
        "cocoa": revision(directory=cocoa),
        "core": revision(directory=core),
        "lsst_y1": revision(directory=project),
        "source_nz_sha256": sha256(filename=source_file),
        "lens_nz_sha256": sha256(filename=lens_file),
        "files": {},
    }
    for filename in sorted(output.iterdir()):
        manifest["files"][filename.name] = sha256(filename=filename)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Saved source {args.source_bin}, lenses {args.lens_bins} to {output}")
    print(f"CAMB sigma8 = {sigma8:.12g}; no covariance matrix was computed.")


if __name__ == "__main__":
    main()
