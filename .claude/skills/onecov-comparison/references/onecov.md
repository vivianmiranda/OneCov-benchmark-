# Local OneCovariance source study

Studied 2026-10-05 at OneCovariance commit `311c2cf`. Paths below are
relative to that checkout. These are implementation findings, not measured
agreement with CoCoA. Recheck contracts if the checkout changes.

## Where the work happens

| File | Responsibility |
| --- | --- |
| `covariance.py` | Native INI entry point; constructs the estimator and writes output. |
| `onecov/cov_input.py` | Configuration and input-table parsers; supplies many defaults. |
| `onecov/cov_setup.py` | Cosmology, input consistency, survey and redshift setup. |
| `onecov/cov_halo_model.py` | Mass function through `hmf`, profiles and halo integrals. |
| `onecov/cov_hod.py` | Central/satellite occupation and mass-observable relations. |
| `onecov/cov_polyspectra.py` | Matter/tracer power, responses and halo trispectra. |
| `onecov/cov_ell_space.py` | Radial projections and Gaussian, SSC and cNG angular covariance. |
| `onecov/cov_theta_space.py` | Flat-sky Bessel transforms into angular correlations. |
| `onecov/cov_output.py` | Component lists, combined matrices and saved spectra. |
| `LevinBessel/` | Bundled C++ oscillatory integrator, imported as `levin`. |

COSEBIs, bandpowers and arbitrary summaries have separate estimator
classes. OneCov's photometric/spectroscopic 6x2pt is not Cocoa's cluster
6x2pt plus counts. Its stellar mass function is also a different observable.

## Input contracts that can change the physics

**Source naming.** `zlens_file` means the source population used for
lensing; `zclust_file` means galaxy lenses. Each text file is `z, n1, n2, …`.
`value_loc_in_lensbin` and `value_loc_in_clustbin` distinguish midpoints
from histogram edges. LSST's current files contain midpoints. OneCov may
drop an initial near-zero midpoint; inspect the parsed support and do not
silently alter the input to conceal this difference.

**Noise.** `n_eff_*` is in arcmin^-2. `ellipticity_dispersion` is the
dispersion of one shape component, so use 0.26, not sqrt(2) times 0.26.
Each selected LSST source bin keeps density 2; each selected lens bin
keeps density 3.6. Selecting one bin does not move the entire catalogue
into that bin. These equal per-bin densities belong to Cocoa's forecast,
not a reconstruction of its supplied likelihood covariance.

**Linear bias.** Supply `[bias] bias_files` with columns `z, b1, …`, and
set `[observables] unbiased_clustering = False`. `PolySpectra.__init__`
then sets its internal matter-tracer mode and `redshift_dep_bias = True`.
`CovELLSpace.__set_redshift_distribution_splines` multiplies each galaxy
window by that bin's bias. The numerical value of the internal
`unbiased_clustering` flag alone therefore does not describe the final
projected galaxy bias. A missing bias file can lead to a warning/fallback;
verify parsed tables before a run.

**Amplitude.** OneCov accepts `sigma8`; Cocoa's LSST forecast specifies
`As = 2.1e-9`. Compute sigma8 with CAMB at the same cosmology rather than
guessing it. Matching sigma8 alone does not make internally generated
power tables identical: solver resolution, nonlinear prescription and
interpolation remain independent choices.

**Power tables.** `Pmm_file` uses `z, k, P` with k in h/Mpc and P in
(Mpc/h)^3. All k values for the first redshift precede the second
redshift; use increasing z and at least two distinct z values. These
tables feed projected power. They do not replace every internal linear
power calculation used by the halo response/trispectrum.

**Angular tables.** Prefer the explicit four-column form
`ell, tomo_i, tomo_j, C_ij`. Use 1-based tomography; include all ordered
pairs, including both (i,j) and (j,i), with ell outermost. For one source
bin this is `ell, 1, 1, C`. C is dimensionless and excludes shot/shape
noise. Avoid the two-column shortcut: its recursive branch is not the
same as the documented explicit-table path.

**Angular bins.** `ell_min/max/bins/type` selects the internal spectrum
grid. Separate `ell_*_lensing` and `ell_*_clustering` settings select
output bands; ggl shares clustering bands. Logarithmic band edges are
truncated to unique integers, then centers are geometric means. Match
the actual edges and weighting, not just a quoted number of bins.

**Components.** Set `split_gauss = True` and save a list as well as a
matrix. The list preserves sample variance, mixed signal/noise, pure
noise, cNG and SSC. A matrix file contains their sum. With split disabled,
SSC may be merged into cNG and its own column returned as zero.

## Numerical work and expected cost

The normal entry point imports all estimator classes even for a small
Gaussian run. `hmf`, CAMB, healpy and the bundled `levin` extension are
runtime dependencies. `setup.py` builds `levin`, not an installed Python
package named `onecov`; run from the source checkout or add it to the
Python path. Do not assume an unrelated package called `levin` exposes
the bundled `Levin` API.

`CovELLSpace` inherits halo setup. Shared C_ell tables do not eliminate
that initialization. Its `get_Cells` accepts tabulated spectra, but the
constructor still prepares background/halo objects and CAMB.

Important controls are separate: `delta_z` for power evolution,
`tri_delta_z` for trispectra, `integration_steps` for radial work,
`log10k_bins` under power/trispectrum evaluation, and `M_bins` under halo
evaluation. Initial inexpensive settings are pilots, not accuracy defaults.
Changing a node count need not preserve old nodes; design convergence
checks accordingly rather than assume monotonic error.

`[misc] num_cores` reaches some C++ integrators. Contrary to an old
docstring, `__trispectra_234h` at this revision loops serially over k pairs
using SciPy `quad`; merely requesting eight cores does not parallelize
this loop. It caches angular integrals with growth-factor rescaling.
Measure cost before increasing the trispectrum k grid: the independent
pairs grow as N(N+1)/2. Integration-warning fallback samples 3,000 angles;
retain warnings in the log and never classify such a run as converged
without further checks.

## Physical differences to keep visible

- NLA enters the OneCov projection windows and thus G, SSC and cNG.
  Cocoa currently applies its requested IA model only to G. Initial
  comparisons therefore use zero IA; TATT is not an established shared
  model here.
- `calc_Cells_nonLimber` has gg, gm and mm paths. At this revision its
  actual low-ell cutoffs are 100 for gg/gm and 50 for mm, despite different
  prose in the docstring. It integrates geometric-mean unequal-time power
  with Levin; this is not Cocoa's separable linear FFTLog correction plus
  a nonlinear Limber remainder. Keep the initial comparison Limber.
- Real-space transforms use flat-sky J0/J2/J4 kernels. Cocoa uses full-sky
  bin-averaged kernels. This difference persists even if spectra agree.
- No supplied footprint invokes an analytic circular/top-hat model.
  Check its normalization and long-mode approximation against Cocoa's
  spherical-cap variance before an SSC comparison at 12,300 deg^2.
- Native halo choices include `hmf` mass-function/bias models, a mass
  definition, Duffy concentration normalization and optional one-halo
  damping. The example's defaults are not automatically Cocoa's choices.
- `trispectra` exposes tracer combinations and a one-halo contribution;
  `__trispectra_234h` computes 2h, 3h and 4h locally but returns their sum.
  Separate exports need an explicitly reviewed adapter or upstream API
  change. Do not claim that public output already separates them.

The source cites Pielorz et al. (2010) for higher-halo trispectra; the
OneCovariance methods paper is arXiv:2410.06962. A future equation-level
comparison should check these alongside the Krause/Takada sources used by
CosmoLike. This source study alone is not a paper-level validation.

## Next scientific gates

1. Test parser round trips and one-source shared-spectrum Gaussian output.
2. Check exact band definitions, noise and covariance normalization with
   an independent analytic Gaussian calculation.
3. Test native spectra with the same cosmology, then supplied CAMB power.
4. Match halo choices and inspect single-redshift response/trispectrum
   inputs before drawing conclusions from a projected cNG difference.
5. Run one-bin G+SSC and G+cNG pilots separately, refine one control at a
   time, then add a lens/source cross block. Only then budget full LSST.
