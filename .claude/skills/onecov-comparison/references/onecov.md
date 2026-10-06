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

The one-lens case was tested and fails in the unmodified constructor:
the bias reader returns a 1D array but `CovELLSpace.__init__` indexes it
as `[tomo, z]`. The small galaxy pilot therefore selects **two distinct
LSST lens bins**. Do not duplicate a population or patch the array shape
silently. A two-lens, one-source Gaussian run completed with the native CLI.
The parser also requires an HOD configuration even for shear-only and
supplied-bias runs; the pilot carries the upstream example's ancillary HOD
values. They are not an LSST HOD fit.

**Amplitude.** OneCov accepts `sigma8`; Cocoa's LSST forecast specifies
`As = 2.1e-9`. Compute sigma8 with CAMB at the same cosmology rather than
guessing it. Matching sigma8 alone does not make internally generated
power tables identical: solver resolution, nonlinear prescription and
interpolation remain independent choices.

**Power tables.** `Pmm_file` uses `z, k, P` with k in h/Mpc and P in
(Mpc/h)^3. All k values for the first redshift precede the second
redshift; use increasing z and at least two distinct z values. Retain the
first table node above the maximum n(z) redshift, even if its value lies
outside the physical survey: it supplies the last interpolation bracket.
OneCov checks the full n(z) file support, including zero-density tails. These
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
In this revision `__bin_cov_ell_gauss` averages integer multipoles with
uniform weight. For a single bin [L,U) it sums the Gaussian numerator
divided by (2ell+1), then divides by fsky*(U-L)^2. Cocoa's Fourier operator
uses mode-count weights. Do not compare their native band matrices as
identical estimators. `compare_gaussian.py` supplies a benchmark-only
uniform-weight operator to Cocoa's production Gaussian kernel; neither
code's numerical source is changed.

OneCov linearly interpolates the **Gaussian numerator**, after forming
products of spectra. Interpolating each spectrum before forming products
is a different approximation on a coarse grid. The assembly comparison
therefore exports every integer ell in short intervals (30–150 and
1500–1620 inclusive), with five bands that exclude the final endpoint.
`get_Cells` replaces the internal ell grid with the supplied spectrum grid.
Changing the template's `ell_bins` cannot refine a supplied spectrum.

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
- With no supplied footprint, `Setup.calc_a_lm` builds a **spherical cap**
  using HEALPix `query_disc` at NSIDE=1024, then `anafast(use_weights=True)`.
  The radius is arccos(1-area/(2*pi)); retained multipoles are 0..3070.
  This corrects the earlier description of an analytic circular fallback.
  Cocoa uses analytic spherical-cap harmonics. Check pixel-area versus
  nominal-area normalization and harmonic truncation before comparing SSC.
  In particular, Cocoa checks raw C_0=area^2/(4*pi) to 1e-8: never pass a
  pixelized mask with a mismatched nominal area into this fatal C guard.
  OneCov's `survey_variance_mmmm` omits chi^-2, which appears later in
  its projection; Cocoa's sigma_b^2 includes it. Compare complete formulas,
  not intermediate arrays with different distance factors.
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

1. Done: exported project inputs and tested the native parser/CLI for a
   one-source Gaussian matrix and a two-lens, one-source Gaussian matrix.
2. Done for the single-source pilot: the independent uniform-ell Gaussian
   sum agrees within 3.4e-7, limited by native text output precision.
   This checks assembly, not spectrum interpolation convergence.
   The subsequent shared-integer-ell comparison also calls Cocoa's C kernel
   and checks all 3x2 cross blocks. Both low-ell shear and 3x2, and high-ell
   3x2, reproduce NumPy to OneCov's saved precision for every component.
   Cocoa agrees with NumPy to floating-point precision. All totals are
   positive definite; maximum generalized variance changes are 1.69e-7,
   4.65e-6 and 4.62e-7. The component list has only five significant digits;
   the total matrix has seven. Do not misidentify output rounding as a
   physical discrepancy or confuse this with native spectrum convergence.
   See `results/gaussian_assembly_20261005.json` for the reviewed record.
   `time_gaussian.py` subsequently checks native in-memory components:
   both codes agree with the independent reference to roundoff. The timer
   excludes setup, spectrum generation, writing and OneCov block repacking.
   Cocoa uses three production calls (total, signal, noise), deriving the
   mixed component by subtraction. OneCov uses public `covELL_gaussian`
   with split output; its 18 blocks are reordered only for validation.
   No numerical source is patched. See the separate timing record for
   first-call times, raw batch samples and the variable tiny shear case.
   Never describe the approximately 40x 3x2 assembly ratio as a full-survey
   or end-to-end CLI speed ratio.
3. Native and supplied-CAMB-power paths ran successfully. Their relative
   spectra differences and interpolation convergence still need assessment.
4. Match halo choices and inspect single-redshift response/trispectrum
   inputs before drawing conclusions from a projected cNG difference.
5. One-bin G+SSC and G+cNG pilots ran separately and passed output checks;
   their positive totals do not establish component accuracy. Refine one
   control at a time, then add non-Gaussian lens/source cross blocks.
   Only then budget full LSST.

## SSC response comparison started

`ssc_response.py` exports OneCov's native matter response, P_linear,
I11, I02, I12 and its differentiated spectrum at z=0, 0.5 and 1. Export
runs in the OneCov environment, comparison in Cocoa's environment.
The native one-halo helper supplies I02; `I_alpha_xy(alpha=0)` is not an
implemented I02 path in this revision. No source or methods are patched.

With shared moments, the same P_linear slope and fractional=False, Cocoa's
production response routine agrees with OneCov to 3.91e-16 fractionally.
NumPy checks the expression independently. This validates assembly only,
not independent halo moments or the angular SSC matrix. Reviewed results
are in `results/ssc_response_20261005.json`.

The model choices still differ. OneCov differentiates log(k^3 P_linear).
Cocoa's survey default differentiates log(I11^2 P_linear), using 47/21
instead of 68/21 because its input slope excludes k^3. That is the corrected
Takada-Hu (2013) Eq. 44. Cocoa additionally transfers the fractional halo
response to a nonlinear target power. On shared OneCov ingredients, adding
only the I11^2 slope changes responses by at most about 0.23%, 0.18% and
0.14% at these redshifts over 0.001<=k<=10 h/Mpc. This pilot diagnostic is
not a converged default-to-default comparison and excludes target transfer.

Next: match the survey-window convention and projected shear response,
then assess response/mass-grid and radial convergence separately. Native
matter integrals temporarily extend the lower halo mass to 10^2 Msun/h;
their true integration range cannot be inferred from the INI's M_min alone.

## Shared SSC projection, 2026-10-05

`compare_ssc.py` runs the unchanged OneCov `covELL_ssc`, then feeds its
actual shell response, linear long-mode power, lensing window and Simpson
radial rule to Cocoa's production kernels. `capture_native.py` passively
reads named locals at function return using Python's profiler. It does
not replace methods or alter arithmetic. Its instrumented elapsed times
are diagnostic costs, never an uninstrumented performance comparison.

Six 100x100 shear matrices at ell=30..3000 agree within 6.47e-15 in
|delta C_ij|/sqrt(C_ii C_jj). All are symmetric and positive definite.
Radial grids 300,601,1201 were checked; at delta_z=0.05 the 1201->2401
change is 4.01e-5. At 2401 nodes delta_z=0.05->0.025 changes the matrix
by 4.27e-4. This does not establish halo/power-grid or Fisher convergence.
The native input reader can enlarge an underspecified ell grid; always
archive actual model.ellrange, not merely the requested ell_bins.

OneCov's pixelized cap has monopole area 0.048998% below the nominal
12300 deg^2 area. The comparison supplies that monopole area to Cocoa's
raw-mask guard, then explicitly restores OneCov's nominal area squared
normalization. This adapter is only unit/normalization matching. Replacing
the mask with Cocoa's analytic cap at the same nominal area and L_max=3070
changes the shared SSC by at most 0.0395%. No numerical source is edited.

Results: `results/ssc_projection_20261005.json`; plot and reproduction in
README. The source inputs/results are preserved under work/ssc_projection_*
and work/ssc_comparison_*. Earlier smoke/failure folders remain preserved.

## Native halo comparison

`compare_halo.py` evaluates the native models at z=0.1,0.5,1, avoiding
the core halo reader's excluded a=1 endpoint. OneCov mass nodes were
200,400,800; Cocoa mass GL rules were 96 and 256 per panel. Over
k=0.001..10, the final OneCov refinement changes I11/I02/I12/response by
<=1.17e-5, the Cocoa refinement by <=3.98e-7. Native differences persist.

The bias formula at equal nu agrees to 5.56e-16. OneCov divides its raw
Tinker bias by norm_bias=0.772506,0.711928,0.641533 from a finite range.
Cocoa leaves the raw bias and normalizes f(nu) through int b f dnu=1.
OneCov's I11 multiplies the normalized bias by norm_bias again and adds
an unresolved-mass completion; I12/I13 retain the normalized bias. Do not
assume that matching the label Tinker10 matches the moments.

OneCov uses Duffy08; Cocoa Bhattacharya13. Also record the small density
constant difference: 8.326098817e10 versus 8.325600500e10 in Msun/h per
(Mpc/h)^3. With supplied equal concentration, max absolute NFW difference
is 1.39e-5 over the export; it is not a native-concentration agreement.
Native dP/ddelta_b differs by up to about 15%,21%,26% over the sampled k
range, including concentration, bias normalization and fractional transfer.
No model is silently retuned. Next compare separated trispectrum terms,
using shared inputs to distinguish assembly from these halo choices.

## Bias normalization and separated cNG, 2026-10-06

`diagnose_bias.py` integrates the actual hmf and Cocoa fitted functions
over ln(nu)=-90..3.5, at z=0.1,0.5,1. It is not an independent covariance
implementation. Doubling 4097 to 8193 nodes leaves integrals unchanged to
roundoff. hmf mass normalization gives int f=1 but int b_raw f equals
0.993555,0.974448,0.953907. OneCov's finite M=1e2..1e17 interval instead
has norm_bias=0.772506,0.711928,0.641533. The lower nu limits are
0.233868,0.287245,0.362994, so the cutoff excludes a substantial fraction
of the extrapolated fit. Dividing every bias by the finite integral is
an extra prescription, not the calibrated Tinker Eq.7 constraint itself.

Cocoa retains b_raw, normalizes f for int b f=1 (measured within 1.1e-9),
but then int f=1.006486,1.026222,1.048320. Do not present this as satisfying
both constraints or proof of simulation accuracy. OneCov I11 cancels
norm_bias and adds missing low-mass weight; I12 and I13 retain it.
Thus, at fixed other ingredients, the latter rescaling multiplies
2h(13),3h by 1/N and 2h(22) by 1/N^2. Neither 1h nor 4h changes.

`compare_trispectrum.py` captures unchanged native 234h locals and the
native damped 1h. Three z, nine k=0.001..10, all 45 triangular pairs.
OneCov angular int_2h/int_3h must be divided by two; restore their growth
factors first because the object stores D^-2,D^-4,D^-6 factorizations.
Cocoa's production algebra then agrees for 1h,3h,4h to roundoff. For 1h
this is only a supplied-term copy, not a halo integral validation.

The 2h off-diagonal discrepancy is 11.82%,12.76%,13.64% by redshift.
OneCov's integral_mmm[i,j] is I13(K,Q,Q), but both 13 partitions use it;
the second requires its transpose I13(Q,K,K). Cocoa uses both. Feeding
the repeated moment into Cocoa removes the difference to <3.4e-16.
This diagnostic alters supplied arrays only, never either source.
Takada-Hu Eq.29 gives the four permutations and fixes this bookkeeping.
Do not silently repair OneCov or call the native trispectra matched.

Unequal-pair tree angular averages with the same quadratic log-P reader
agree within 2.8e-7. Equal pairs are excluded from that particular test:
an initial diagnostic extrapolated OneCov's quadratic log-P below its
1e-5 table limit and grew unphysically at k approaching zero. The failed
shared-domain experiment is preserved under trispectrum_cocoa_96_v2 and
trispectrum_cocoa_800_*. The final script guards the internal k domain and
compares only unequal pairs; native predictions still include diagonals.
Native corner controls 1e-3->1e-5->1e-7 change some diagonal 2h/3h entries
by tens of percent. Do not claim angular convergence there. Native Cocoa
96->256 mass/angular rules change all sampled terms by <3.93e-6.

Native model curves intentionally retain Duffy/Bhattacharya, bias and
damping differences. All quantities passed to algebra/projection kernels
use consistent Mpc/h units; native Cocoa output converts by (2997.92458)^9.

`compare_connected.py` passively captures native covELL_non_gaussian while
public calc_covELL constructs five-band G+cNG. Captured unbinned cNG has
no survey-area divisor yet; the later native binning adds 1/area. Cocoa
receives the actual trispectrum at eight multipoles, native W^2 and
Simpson dchi/(area*chi^6), in one generic transform slot. Compare all 64
entries. The native 5-band total separately supplies positivity and
generalized variance-ratio diagnostics; never conflate these matrix sizes.

Source code remains unchanged. All measured runtimes in these captures
include instrumentation and are diagnostic costs, not performance claims.

Ten cNG runs: k9,k17,k33,k65,k129 at radial300,dz0.5,mass400;
radial601 at k65; dz0.25 then0.125; mass800 and corner1e-5 separately
at k65,radial601,dz0.25. Default corner1e-3, k range1e-4..1e2.
Shared projection agrees within 2.30e-15 variance-scaled in every run.
All native five-band G+cNG totals are positive definite (minimum
correlation eigenvalue >0.957). Final k65->129 change: 2.97% cNG,
0.0903% total modes. Final dz0.25->0.125: 3.15% cNG, 0.311% total modes.
Mass400->800: 0.00130% cNG, 0.000128% modes. Corner1e-3->1e-5:
2.95% cNG, 0.0884% modes. Component changes use the finer cNG diagonal
rms product. Generalized ratios use the finer total as denominator.

`collect_ng_results.py` aggregates the preserved work folders into the
three dated result records; `plot_connected.py` makes scientific figures.
No timing bars or independent reference-library comparisons in README.
Next gates: settle the 2h partition and bias-prescription interpretation
before native equality claims; refine redshift interpolation further and
validate near-diagonal angular treatment/domain independently. Do not
certify full LSST or Fisher convergence from these one-source tests.
