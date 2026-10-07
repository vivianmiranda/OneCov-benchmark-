"""Time matched halo and non-Gaussian projection stages, sequentially.

Projection inputs are the validated exports, already resident in memory.
OneCov's integration statements are extracted unchanged from its source;
the enclosing generation of halo/response tables is outside this timer.
Halo timing instead calls the native methods of both codes at the same
three redshifts and 200 k samples, retaining their different prescriptions.
These are stage timings, never full-covariance or full-survey speed ratios.
"""

import argparse
import ast
import configparser
import contextlib
import hashlib
import importlib.metadata
import json
import os
import platform
import signal
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.integrate import simpson

from common import power_table_record, revision, sha256


def native_projection(onecov, stage):
    """Compile the actual native integrand and Simpson call unchanged.

    The two statements live inside long methods rather than a callable
    projection API. Read their syntax trees from that checkout, preserving
    numerical expressions. Only replace the destination-array assignment
    by a return: allocating/repacking the tomography container is excluded.
    There is no profiler, source patch, method replacement or rewritten
    covariance formula in a timed call.
    """
    filename = onecov/"onecov/cov_ell_space.py"
    source = filename.read_text()
    tree = ast.parse(source)
    method = "covELL_ssc" if stage == "ssc" else "covELL_non_gaussian"
    output = "SSCELLmmmm" if stage == "ssc" else "nongaussELLmmmm"
    function = next(node for node in ast.walk(tree)
                    if isinstance(node, ast.FunctionDef) and node.name == method)
    assignment = next(node for node in ast.walk(function)
                      if isinstance(node, ast.Assign)
                      and isinstance(node.targets[0], ast.Subscript)
                      and ast.unparse(node.targets[0].value) == output
                      and isinstance(node.value, ast.Call)
                      and ast.unparse(node.value.func) == "simpson")
    statements = []
    if stage == "ssc":
        statements.append(next(node for node in ast.walk(function)
                               if isinstance(node, ast.Assign)
                               and ast.unparse(node.targets[0]) == "integrand"
                               and "Pmm_response" in ast.unparse(node.value)))
    statements.append(ast.Return(value=assignment.value))
    body = ast.parse("def project():\n    pass\n")
    body.body[0].body = statements
    ast.fix_missing_locations(body)
    namespace = {"np": np, "simpson": simpson}
    exec(compile(body, str(filename), "exec"), namespace)
    expressions = "\n".join(ast.unparse(node) for node in statements)
    provenance = dict(source_sha256=sha256(filename), expressions=expressions,
                      expressions_sha256=hashlib.sha256(expressions.encode()).hexdigest(),
                      onecov=revision(onecov))
    return namespace, provenance


def projection_call(args, data):
    """Prepare common radial geometry; timed calls allocate fresh results."""
    chi, window = data["chi"], data["window"]
    if args.stage == "ssc":
        reference = data["onecov"]
    else:
        # Native covELL_non_gaussian returns this integral before the
        # later band-averaging stage divides by area. Time that same scope.
        reference = data["onecov"]*float(data["area_sr"])

    if args.backend == "onecov":
        namespace, provenance = native_projection(args.onecov.resolve(), args.stage)
        namespace.update(self=SimpleNamespace(los_integration_chi=chi),
                         weight=window**4/chi**6, i_sample=0, j_sample=0)
        if args.stage == "ssc":
            namespace.update(Pmm_response=data["response"],
                             survey_variance_mmmm=data["variance"][None, None, :])
        else:
            namespace["trispec_integrand_mmmm"] = data["trispectrum"][:, :, :, None, None]
        return namespace["project"], reference, provenance

    import cosmolike_lsst_y1_interface as ci
    backend = ci.covariance
    provenance = dict(interface_sha256=sha256(ci.__file__))
    if args.stage == "ssc":
        response = np.ascontiguousarray(data["response"].T)
        pairs = np.ascontiguousarray(np.broadcast_to(window**2, response.shape))
        mean = np.zeros_like(pairs)
        signal = np.zeros(len(data["ell"]))
        weights = data["dchi"]*data["variance"]/chi**2

        def compute():
            # The shell routine adds the lensing windows and inverse chi^2.
            # Keep it inside the timer; response values are the shared input.
            shell = backend.covariance_ssc_shell_response(
                distance=chi, signal=signal, pair_window=pairs,
                mean_window=mean, power_response=response)
            return backend.covariance_project(left=shell, right=shell,
                                               weight=weights)
    else:
        size = len(data["ell"])
        projected = np.zeros((4*size, 4*size, len(chi)))
        projected[:size, :size] = data["trispectrum"].transpose(1, 2, 0)
        probes = np.array([0], dtype=np.int32)
        pairs = np.ascontiguousarray(window[None, :]**2)
        weights = data["dchi"]/chi**6

        def compute():
            return backend.covariance_project_connected(
                probes=probes, pair_window=pairs,
                projected=projected, measure=weights)
    return compute, reference, provenance


def halo_call(args, data):
    """Request native I11, I12, both I13 partitions and undamped I04.

    All requested k pairs are included. Cocoa additionally computes I02
    in its combined API; that work remains in the timer. OneCov returns
    full square pair arrays, Cocoa their symmetry-compressed equivalents.
    Initialization and resident linear-power/variance inputs are excluded.
    """
    wave = data["k"][0]
    if not np.all(data["k"] == wave):
        raise ValueError("halo benchmark needs a common k grid at all redshifts")
    if args.backend == "cocoa":
        root = Path(__file__).resolve().parents[2]/"cocoa/Cocoa"
        core = root/"external_modules/code/cosmolike_core"
        sys.path[:0] = [str(core), str(root/"projects/lsst_y1/covariance")]
        import cosmolike_lsst_y1_interface as ci
        from lsst_y1_covariance import configuration, initialize
        settings = configuration(gaussian={"nonlimber": False, "ia": "none"},
                                 integration_accuracy=args.integration_accuracy)
        tables = initialize(ci, settings)
        scale = np.ascontiguousarray(1/(1+data["redshift"]))
        wavenumbers = np.ascontiguousarray(data["k"]*2997.92458)

        def compute():
            return ci.covariance.covariance_halo_moments(
                a=scale, k=wavenumbers, lnm_edges=settings["lnm_edges"],
                nquad=settings["halo_mass_nquad"])
        return compute, dict(core=revision(core), interface_sha256=sha256(ci.__file__),
                             installed_power=power_table_record(tables, settings),
                             mass_nquad=settings["halo_mass_nquad"],
                             mass_panels=len(settings["lnm_edges"])-1)

    sys.path.insert(0, str(args.onecov.resolve()))
    from onecov.cov_input import FileInput, Input
    from onecov.cov_polyspectra import PolySpectra
    config = configparser.ConfigParser()
    config.read(args.inputs/"onecov.ini")
    # The native reader ignores the trispectrum settings unless enabled.
    # We call only the moment helpers, not the complete angular trispectrum.
    config["covariance terms"]["nongauss"] = "True"
    config["output settings"]["directory"] = str(args.output.resolve())
    config["trispec evaluation"].update(
        log10k_bins=str(len(wave)), log10k_min=str(np.log10(wave[0])),
        log10k_max=str(np.log10(wave[-1])), small_k_damping_for1h="none")
    config_file = args.output/"onecov.ini"
    with config_file.open("w") as stream:
        config.write(stream)
    terms, obs, output, cosmo, bias, ia, hod, survey, prec = (
        Input().read_input(config_name=str(config_file)))
    tables = FileInput(bias).read_input(config_name=str(config_file))
    models = [PolySpectra(float(z), terms, obs, cosmo, bias, hod, survey,
                          prec, tables) for z in data["redshift"]]
    for model in models:
        if not np.allclose(model.krange_tri, wave, rtol=1e-12, atol=0):
            raise ValueError("native 1h and pair-moment k grids differ")

    def compute():
        result = []
        for model in models:
            i11 = model.halo_model_integral_I_alpha_x(bias, hod, prec["hm"], 1, "m")
            i12 = model.halo_model_integral_I_alpha_xy(bias, hod, prec["hm"], 1, "m", "m")
            i13 = model.halo_model_integral_I_alpha_mmm(bias, hod, prec["hm"], 1)
            i04 = model._PolySpectra__trispectra_1h(
                bias, hod, prec["hm"], False, False, False, False, False, True)[5]
            result.append((i11, i12, i13, i04))
        return result
    return compute, dict(onecov=revision(args.onecov), mass_nodes=prec["hm"]["M_bins"],
                         config_sha256=sha256(config_file),
                         halo_source_sha256=sha256(args.onecov/"onecov/cov_halo_model.py"),
                         polyspectra_source_sha256=sha256(args.onecov/"onecov/cov_polyspectra.py"))


def halo_arrays(result, backend, data):
    """Convert outputs for checks only; no rearrangement enters the timer."""
    first, second = np.triu_indices(data["k"].shape[1])
    if backend == "cocoa":
        single, moments = result
        return dict(i11=single, i12=moments[1]*2997.92458**3,
                    i13=moments[2]*2997.92458**6,
                    i13_exchanged=moments[3]*2997.92458**6,
                    i04=moments[4]*2997.92458**9)
    return dict(
        i11=np.array([r[0][:, 0] for r in result]),
        i12=np.array([r[1][first, second, 0, 0] for r in result]),
        i13=np.array([r[2][first, second, 0] for r in result]),
        i13_exchanged=np.array([r[2][second, first, 0] for r in result]),
        i04=np.array([r[3][first, second, 0, 0] for r in result]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("ssc", "connected", "halo"))
    parser.add_argument("inputs", type=Path)
    parser.add_argument("--backend", choices=("onecov", "cocoa"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--onecov", type=Path, default=Path(__file__).resolve().parents[2]/"OneCovariance")
    parser.add_argument("--repeats", type=int, default=11)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--integration-accuracy", type=int, default=2,
                        help="Cocoa halo mass quadrature level (default: 2)")
    parser.add_argument("--halo-reference", type=Path,
                        help="previous Cocoa halo comparison directory")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.output.exists() or args.repeats < 2 or args.batch_size < 1:
        parser.error("choose a new output, >=2 repeats and positive batch size")
    if int(os.environ.get("OMP_NUM_THREADS", "0")) < 1:
        parser.error("set OMP_NUM_THREADS before starting")
    if os.environ.get("OPENBLAS_NUM_THREADS") != "1":
        parser.error("set OPENBLAS_NUM_THREADS=1 before starting")
    args.output.mkdir(parents=True)
    signal.alarm(args.timeout)
    metadata = json.loads((args.inputs/"export.json").read_text())
    if sha256(args.inputs/"inputs.npz") != metadata["inputs_sha256"]:
        raise ValueError("input export changed")
    with np.load(args.inputs/"inputs.npz") as source:
        data = {key: source[key] for key in source.files}

    # Initialization, source extraction and disk reads finish before timing.
    # The first evaluation is recorded separately; subsequent calls still
    # recompute outputs. None of these numbers is a full-construction time.
    with (args.output/"initialization.log").open("w") as log:
        with contextlib.redirect_stdout(log):
            started = time.perf_counter()
            if args.stage == "halo":
                compute, provenance = halo_call(args, data)
            else:
                compute, reference, provenance = projection_call(args, data)
            setup = time.perf_counter()-started
            started = time.perf_counter()
            first_result = compute()
            first_seconds = time.perf_counter()-started
    samples = []
    for repeat in range(args.repeats):
        started = time.perf_counter()
        for _ in range(args.batch_size):
            result = compute()
            del result
        samples.append((time.perf_counter()-started)/args.batch_size)
    final_result = compute()

    if args.stage == "halo":
        first = halo_arrays(first_result, args.backend, data)
        final = halo_arrays(final_result, args.backend, data)
        checks = {name: float(np.max(np.abs(final[name]/first[name]-1)))
                  for name in first}
        if not all(np.isfinite(value).all() for value in final.values()):
            raise ValueError("nonfinite halo result")
        if max(checks.values()) > 1e-12:
            raise ValueError("halo inputs changed during repeated calls")
        # The earlier halo study independently called these same native APIs.
        # Verify the diagonal subset after timing, using its archived results.
        if args.backend == "onecov" or args.halo_reference is not None:
            baseline = data if args.backend == "onecov" else np.load(
                args.halo_reference/"cocoa.npz")
            left, right = np.triu_indices(data["k"].shape[1])
            checks["saved_i11_max_fractional"] = float(np.max(np.abs(
                final["i11"]/baseline["i11"]-1)))
            checks["saved_i12_max_fractional"] = float(np.max(np.abs(
                final["i12"][:, left == right]/baseline["i12"]-1)))
            if max(checks.values()) > 1e-11:
                raise ValueError("halo moments differ from saved native results")
        np.savez_compressed(args.output/"outputs.npz", **final)
        shape = [len(data["redshift"]), data["k"].shape[1]]
    else:
        rms = np.sqrt(reference.diagonal())
        checks = dict(max_variance_scaled_residual=float(np.max(np.abs(
            final_result-reference)/rms[:, None]/rms[None, :])),
            repeated_max_absolute=float(np.max(np.abs(final_result-first_result))))
        if not np.isfinite(final_result).all() or checks["max_variance_scaled_residual"] > 1e-11:
            raise ValueError("timed projection disagrees with native saved matrix")
        np.savez_compressed(args.output/"outputs.npz", matrix=final_result)
        shape = list(reference.shape)
    record = dict(stage=args.stage, backend=args.backend, shape=shape,
                  radial_nodes=len(data["chi"]) if args.stage != "halo" else None,
                  mean_seconds=float(np.mean(samples)),
                  std_seconds=float(np.std(samples, ddof=1)),
                  median_seconds=float(np.median(samples)), samples_seconds=samples,
                  first_call_seconds=first_seconds, setup_seconds=setup,
                  batch_size=args.batch_size, checks=checks, implementation=provenance,
                  input_export=metadata, script_sha256=sha256(__file__),
                  output_sha256=sha256(args.output/"outputs.npz"),
                  versions={name: importlib.metadata.version(name)
                            for name in ("numpy", "scipy")},
                  python=sys.version,
                  platform=platform.platform(), thread_environment={name: os.environ.get(name)
                  for name in ("OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES",
                               "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")})
    (args.output/"timing.json").write_text(json.dumps(record, indent=2)+"\n")
    print(f"{args.backend} {args.stage}: {record['mean_seconds']*1000:.4f} "
          f"± {record['std_seconds']*1000:.4f} ms; checks passed")


if __name__ == "__main__":
    main()
