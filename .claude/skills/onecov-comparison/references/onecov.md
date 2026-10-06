# Local OneCovariance source study

## Complete shear cross-code figure

`complete_shear.py` computes a coherent 100x100 G/SSC/cNG/total pilot for
LSST Y1 source bin 3. `plot_complete_shear.py` plots CoCoA minus OneCov;
all panels use the OneCov total diagonal rms product, with percent units.
The README embeds results/figures/complete_shear_difference.png near its
top. This is not the complete 1560-entry real-space survey calculation.

Use ell=30,60,...,3000 and the band-centre approximation with delta_ell=30.
Removing OneCov's output-band controls makes its public calc_covELL use
the linear-grid Gaussian mode width of 30. CoCoA's production forecast
uses single-integer-multipole operators, then G alone is divided by 30.
SSC/cNG are centre samples and receive no width factor. Neither code's
numerical source or callable is modified. Do not describe this estimator
as an integrated broad-band covariance. It avoids the different native
Gaussian broad-band weights rather than silently comparing them.

Fresh inputs work/lsst_y1_complete_20261006 supply the same CAMB nonlinear
power for Gaussian projections, with each code retaining its native halo
and SSC prescriptions, internal linear-power preparation and footprint.
The known amplitude/radiation, halo and 2h/corner differences below remain.
CoCoA keeps its adopted 1e4 mass minimum; OneCov keeps 1e6. This comparison
is not another shared-response or shared-trispectrum assembly test.

Results: work/complete_shear_{onecov_v2,cocoa}_20261006, summarized in
results/complete_shear_20261006.json. OneCov uses 601 radial samples, 129 k
samples, 400 mass samples, delta_z=.05 and tri_delta_z=.125. CoCoA uses
the LSST default AB1/integration0. All entries are retained, totals are
positive definite. Maximum absolute differences in percent of the total
rms product: G .688563, SSC .241387, cNG .226250, total .328315. Total
variance ratios span .9540080430..1.0069607884. Native numerical/Fisher
convergence is not established by this pilot. Do not conflate these
cross-code differences with error estimates for either implementation.

Single sequential construction timings: OneCov 70.171 s, CoCoA 20.704 s;
setup is separate at 11.673 s and 0.419 s, respectively.
These centre-sample timings cannot predict the full real-space transforms.
The first attempted export used the older lsst_y1_small power table, which
stopped at z=3 and correctly failed OneCov's z-support guard. Its output
folder/log was preserved; fresh prepare_lsst_y1.py output includes the
bracketing CAMB node beyond the source support and resolves that input error.

## Real-space cross-code figure, 2026-10-06

User explicitly requires OneCov's own real-space numerical settings, not
controls chosen to mimic CoCoA's full-sky cutoff. Use the shipped
config_files/config_3x2pt_rcf.ini with only physical survey/IA/calibration,
selected observables, requested components and output locations adapted.
The saved INI and template fingerprint record exactly what was used.
Its native ell=2..1e5 (500 log samples), 500 radial nodes, 900 mass nodes,
100 trispectrum k samples, dz=.08, tri_dz=.5, theta_accuracy=.01 and
integration_intervals=400 were unchanged. The Bessel weight support is
hardcoded 1..1e5 in this revision; distinguish that from spectrum samples.
CCL-benchmark's original 0044cc6 LSST FFTLog script supplied ell<=30000;
never transfer a full-sky cutoff requirement to a different transform.

`complete_shear.py --space real` computes source bin 3 at eight common
logarithmic annuli 2.5..250 arcmin, ordered xi+ then xi-, 16x16. It calls
native CovTHETASpace.calc_covTHETA and the CoCoA production interface.
OneCov split-G groups 7,8,9 (three entries each) and NG groups 7,8,9 are
++, +-, --. Reverse ordering is the transpose of +-. Native auto-block
asymmetries remain untouched. No code in either numerical library changed.

Archives: work/real_shear_onecov_native_20261006 and
work/real_shear_cocoa_20261006. Report results/real_shear_20261006.json;
PNG/PDF results/figures/real_shear_difference. Every entry is included.
Max percent difference / OneCov total rms product: G2.018179,
SSC0.662422, cNG0.267596, total2.124543. Generalized total mode ratios
.9847814351..1.0234774897. Symmetric-part correlation minimum eigenvalues
are .16254508 (CoCoA), .16124295 (OneCov). Native OneCov maximum asymmetry
is 1.13565e-5 of total rms (G); CoCoA is exactly symmetric. Eigenvalues
use explicit symmetric parts only; plots and NPZ preserve raw matrices.
Do not call these repaired or claim exact symmetry for OneCov.

Sequential single runs at8 threads: CoCoA setup .419010s + construction
45.443902s; OneCov setup19.910623s + construction513.642171s. These are
native-setting construction measurements, not equal-accuracy timings.
No additional numerical job ran simultaneously. Documentation and Git
checks ran while waiting. Full LSST OneCov has not been timed. Hours or
longer is a planning estimate, not a validated bound from matrix size.
Native cNG convergence and attribution of the 2% Gaussian difference are
pending; no Fisher-convergence claim. Do not retune OneCov defaults silently.

## Simultaneous halo mass and bias normalization, 2026-10-06

The later request authorizes sensitivity tests, not adoption of a new
production normalization. `halo_normalization.py` supplies an explicit
backend adapter; it does not replace installed functions or change either
code checkout. `run_normalization_covariance.py` uses the native LSST YAML
and production backend. Choices: native, mass_only (f/S, b), mass_and_bias
(f/S, S*b), and ingredients. S is the analytic integral of the actual
CoCoA Tinker shape and interpolated amplitude over all peak heights.

The two Gamma-function shape integrals agree with numerical integrals of
the installed fnu to better than 1e-10 over z=.01,.1,.5,1,2,3. Public fnu
requires a<1; do not test it at a=1. The scalar NFW formula is singular at
k=0; use the analytic unit limit, as the production covariance does.
Direct mass-quadrature checks at z=.1,.5,1 and k=0,.01,1,100 h/Mpc verify
I11 and all five pair roles, including off-diagonals, within 2e-10 relative.
Mass-only I11 must be u_min+(I11-u_min)/S, not I11/S: recompute completion.
For paired normalization only I02/I04 change; n*b stays unchanged.

Full 1560x1560 archives, all entries/no cuts:
work/normalization_native_v2_20261006,
work/normalization_mass_only_20261006,
work/normalization_mass_and_bias_20261006.
Identical power hashes, installed interface, survey grids and settings.
Comparisons and four-component PNG/PDF panels are in results/normalization*.
Max entry differences / native total rms, percent (G, SSC, cNG, total):
mass-only (0,.0361054,.8581820,.8677758);
paired (0,1.1669957,.5519169,.6245842).
Generalized total ratios: .97713295394..1.00007166622 and
.98373668873..1.01607758359, respectively. All totals positive definite;
Gaussian is bitwise unchanged. Runs took 51.07,50.90,51.04 seconds, each
sequential with eight OpenMP threads/one BLAS; timings include adapter
overhead and are not a production speed comparison. No Fisher or improved
physical accuracy claim follows from these sensitivity tests.

The first native archive and the v2 post-construction check failures were
preserved. They concern diagnostic endpoint handling, not matrix generation.
V2 saved the complete matrix and report before the failed zero-profile
check. The corrected ingredient-only run passed separately, without
rerunning or replacing that matrix. Both trial logs remain in /tmp.

Physics anchors: Tinker 2010 Eq.7 and text after Eqs.9--12 explicitly
sets alpha by the biased integral; no missed numerical division is implied.
HMcode-2020 (2009.01858), Sections3.1/3.2/4.4, adopts mass-normalized ST
and replaces the standard two-halo profile integral by its unit limit in
the production matter power. It is not a covariance calibration.
An unresolved component with mass 1-F and bias (1-B)/(1-F) can satisfy
both finite-population constraints while preserving resolved fits. This
interprets the existing completion weight 1-B; it is not an implemented
new population with all higher moments or a validated low-mass model.
Do not promote a global bias rescaling solely because both integrals close.

## Sub-solar cutoff diagnostic, 2026-10-06

Production remains at Mmin1e4. The isolated builder can extend the table
domain to1e-3 and sigma kmax to1e7/1e9, retaining the existing log-log power
continuation. No installed source or library changes. The lower negative
mass edge must use math.log(10**exponent), matching the scalar C guard;
exponent*log(10) can round below its strict domain by one ulp.

Study work folders: microhalo_k5_boost4_v2, microhalo_k7_boost4,
microhalo_k9_boost4 and microhalo_k9_boost8. Failed first k5 folder is
preserved. Four cutoffs1e4/1e2/1/1e-3, z=.1/.5/1,41k=.001..100,861pairs.
Mass quadratures96/128/256/512 were all tested; do not skip128 in future
scans. Accuracy refinements are validation, not automatic demands for
production precision. Native96 differs from512 by<=3.83e-5 fractionally
over the tested moments/terms/response.128 gives5.40e-6;256 gives2.70e-7.
Fine table boost4->8 changes ingredients<=9.46e-6, but the cutoff contrast
changes<=9.33e-11. This does not certify full Fisher/physical accuracy.

At z1, missing bias response1e4->1e-3 drops .28384818->.20007407;
ordinary resolved mass .60534547->.74645495. The full-fit mass integral
S=1.04832 is unchanged:1-F is a deficit relative to unit matter density,
not the literal integral of the fitted unresolved population S-F.

Sigma kmax1e5 underestimates the smallest mass sigma by~5.3% relative to
extended FFTLog;1e7->1e9 changes it<1.91e-9. OneCov's unchanged hmf TopHat
on the same continued P agrees with fine sigma at the four masses within
1.92e-6. About90.14% of sigma^2 at1e-3 comes from P above suppliedkmax143.
Numerically extending the reader is possible; a calibrated microhalo/DM
free-streaming model is a distinct physical question.

Full1560 matrices: work/microhalo_full_native4, microhalo_full_wide4,
microhalo_full_deep. Common-wide4->deep changes total modes1.869e-9,
G bitwise unchanged, both totals positive. Sigma-domain/kmax expansion
alone atMmin1e4 changes modes1.76824e-5. All four complete component
plots and diagnostics are saved as results/microhalo*. Construction
times52.126,51.912,52.362s: single sequential8-thread measurements with
background application activity, not a precise slowdown estimate.

The authorized split-rule and sequence follow-up is completed below.
No production adoption is authorized by these studies.

## Split low-mass quadrature and sequence acceleration, 2026-10-06

`split_mass_rule.py` changes only quadrature construction in a private
copy of halo_cov.c. Upper panels retain the public nquad. Panels below
1e4 use COCOA_DIAGNOSTIC_TAIL_NQUAD, set explicitly by benchmark runners.
Native weights, profiles, SIMD sums and one final completion are retained.
The installed sources and library are untouched.32 is an explicitly
authorized diagnostic exception, not a change to the production >=64 rule.
The private build uses mass minimum1e-20 and sigma kmax1e15, retaining the
same log-log power continuation. Preserve all original upper mass nodes.

Folders: work/split_tail_build, split_single32, split_intervals32,
split_intervals128, split_refined. Upper quadrature96 throughout. Single32
uses one24-decade lower interval. Intervals use six4-decade panels with
32 or128 nodes. Refined uses128 lower nodes and table boost8 rather than4.
Partial sums at exponents4,0,-4,-8,-12,-16,-20 retain every earlier node.
All native corrected zero modes and repeatability checks pass.

Single32 versus intervals128 changes the tested corrected ingredients
<3e-12 fractionally. Missing z1 response at1e-20 is.09207 and resolved
ordinary mass.92024. Fine table changes are~9.46e-6, much larger than the
quadrature contrast; never call the latter absolute physical accuracy.
At the smallest mass, sigma differs from the same-power OneCov hmf filter
by.116%. The filter's1e13->1e15 tail change is1.11e-10.99.84% of variance
comes from power above suppliedkmax143. Extreme-domain FFTLog/table
accuracy and small-halo calibration remain limits of this diagnostic.

Full folders: split_full_base4, split_full_single32, split_full_reference128.
Times49.168,50.702,50.406s, sequential8 threads, single runs with background
app activity: no precise speed claim. All1560 entries in G/SSC/cNG/total
are checked. G is bitwise unchanged and totals positive. Single32 versus
ref128: normalized entry differences<4e-13; mode ratios differ<7e-12 at
eigensolver precision. Tail extension4->-20 on common tables changes total
modes1.8673e-9. Extreme table-domain expansion at fixed1e4 instead gives
8.31318e-5. Component plots/reports are results/split_full_*_20261006.json.

`report_split_tail.py` tests Aitken delta-squared and mpmath's native
Wynn/Shanks (randomized=False).40-digit arithmetic avoids additional
algorithm rounding; inputs remain finite-accuracy double-precision
quadratures. No noise recovery or randomized singularity handling.
Known geometric-series check passes. At z1, raw last partialB=.907929;
Aitken limit1.002344 and Wynn.999010. Wynn with5 partials gives1.000481,
so convergence is not monotone. Held-out next-two partial sums after a
fit through1e-12 agree within.0242%/.0749% acrossz. Refinement128/table8
retains the behavior. All test extrapolations take~.029s in one batch.

Extrapolating raw finite-k I11 after removing the native completion gives
max errors.519% Aitken and.216% Wynn against the deep corrected reference.
Existing completion atMmin1e4 differs only.00729% overk<=100. Thus no
extrapolated forecast replaces that already-better completed moment;
the full matrices test split quadrature, not an accelerated physical model.
For tiny halos u~1, completed I11~1+integral dw*(u-1), explaining why rawB
converges slowly while corrected finite-k moments barely change.

Results, scientific plots and human explanations are in the README and
results/split_tail_20261006.json. Timing comparisons remain tables.

## Deeper acceleration and numerical FFTLog weighting, 2026-10-06

The follow-up goal is to find the shallowest reliable cutoff before
numerical error dominates. Do not dismiss this goal merely because the
existing completion already imposes I11(0)=1. No production adoption or
change to either halo fit is authorized by the diagnostic study.

Private build work/deep_series_build extends Mmin to1e-50, sigma kmax
to1e25, and changes the FFTLog weighting b from1.5 to0.5. The factor
k^-b removed from Delta² is restored by the analytic Mellin kernel.
This changes finite-FFT conditioning/periodic approximation, not physical
power, halo bias or multiplicity normalization. Do not attribute all
improvement specifically to floating-point roundoff without isolating it.

The one-variable control work/split_bias05_control versus split_refined
keeps Mmin1e-20, kmax1e15, table boost8, tail128 and all supplied P fixed.
Only b changes. Native same-power OneCov filter discrepancy at1e-20 falls
from0.116% to0.000173%. report_fftlog_bias.py verifies shared arrays,
domains and controls. Deep same-power sigma at1e-50 differs0.000237%;
the filter kmax1e23->1e25 change is<1.11e-10 fractionally.

work/deep_series32 and deep_series128 use exponents4,0,-4,...,-48,-50;
upper96 nodes unchanged, lower32/table4 versus lower128/table8. The
panel boundaries stay anchored at4,0,-4,...; -50 adds two decades without
moving previous nodes. Do not feed that last unequal interval into the
equal-step sequence. A held-out geometric prediction of it errs<0.000114%.

report_deep_series.py compares Aitken, Wynn last7 and highest available
odd-order window. Highest-order includes the newest sum, dropping only
the oldest when an odd input count is required. Atz1 rawB at1e-48=.974640;
Wynn=.99999615, refined=.99999636. At1e-40, maximum finite-k I11 error
versus fine completed native I11 is0.00668% (fine0.00665%). At1e-48 it is
0.000858% (fine0.000776%). First tested cutoffs passing0.1%,0.01%,0.001%
are1e-32,1e-40,1e-48. This is an ingredient diagnostic at three redshifts
and41 k=.001..100, not a universal covariance accuracy target. Last7 Wynn
has an unstable near-pole at1e-28; more depth/order alone is not a test.

Crucial distinction: the full fitted ordinary mass integral atz1 is
1.04832034, independent of the covariance cutoff. Lowering the cutoff or
accelerating cannot turn it into1. With a unit-mass shape f0, full
integral b*f0=.95390689; Cocoa rescales f=f0/.95390689 to impose unit
bias integral, hence full mass1.04832. This precedes mass truncation.
Never suggest the4.83% mass excess is omitted low-mass weight. The earlier
phrase 'even when the full model satisfies the condition' meant ONLY
the bias-weighted condition. README now explicitly separates both.
At1e-50 directF(z1)=1.019583; Wynn predicts1.0483204. Limits are measured
from the same private build at fine settings in work/deep_mass_limit,
using the analytic mass_integral helper. No normalization is changed.

Full native completed LSST runs deep_full20/deep_full50 take50.297/51.300s
sequentially, OMP8/BLAS1. Every1560x1560 entry compared; G bitwise equal,
SSC/cNG/total changes over reference total rms <=5.45e-14, both positive,
generalized modes <=2.61e-12 from unity (comparison precision). These
test cutoff extension, NOT an accelerated replacement covariance.
No full accelerated forecast has been generated; that validation remains
separate if replacing the existing completion is later requested.

Results and figures: deep_series, deep_series_power, fftlog_bias,
fftlog_bias_control, deep_full_extension, all suffixed20261006.json.
All native source files and the installed library remain unchanged.

Keep the main README's low-mass discussion short: compare original1e6
with the experimental1e-40+Wynn case. Detailed intermediate scans and
reproduction moved to docs/low_mass_study.md. Do not call Wynn the released
or production implementation. The saved original-domain mass_cutoff_native8
control (upper256, table8) has maximum completed-I11 discrepancy0.18046%
from the same deep reference, versus0.00668% for the proposed Wynn case.
Supplied P is identical, but sigma domain/weighting and quadrature also
differ; this is a configuration comparison, not an isolated Wynn speedup.

## Low-mass cutoff study: scope

The requested study reduces the finite-range I11 completion by integrating
more low-mass halos explicitly. Keep CoCoA's multiplicity normalization,
fitted bias, concentration and additive minimum-mass-profile prescription
unchanged. This is not a request to reconcile the two codes by changing
either halo model or to adopt a different normalization convention.

The initial study used an isolated LSST build with only its halo-table
lower mass boundary extended from 1e6 to 1e2 Msun/h, leaving the installed
interface and repository C/C++ sources unchanged. Compare 1e6, 1e4 and 1e2
using the same wide sigma table, preserving every original mass panel.
Refine GSL mass quadrature and internal tables separately. Record missing
bias-weighted response, corrected moments, SSC response, separate 1h--4h
terms, and repeated ingredient times. Check the original-domain build at
the original cutoff to identify mass-table regridding effects.

The supplied CAMB table ends near 143 h/Mpc; the native edge continuation
is held fixed. Check the sigma integration tail using that same continued
power. Extending its numerical range does not physically calibrate the
extrapolated small-scale spectrum or low-mass halo fits. Do not infer
survey/Fisher convergence from the ingredient study or interpret a large
completion weight as an equally large prediction error.

Completed record: `results/mass_cutoff_20261006.json`. At z=0.1,0.5,1,
lowering M_min from 1e6 to 1e2 changes missing weights from
[0.216486,0.262624,0.319728] to [0.166180,0.205728,0.255131]. At z=1 this
is a 20.2% reduction of the correction, not a 20.2% covariance change.
Over k<=10 h/Mpc, corrected I11 changes <=0.00111% and the summed
trispectrum <=0.00123%; over k<=100 these are 0.1795% and 0.2866%.
The 1e4-to-1e2 changes are about 25 times smaller than 1e6-to-1e4.
The isolated 4h term can change 0.7161% over the larger k interval.
Corrected k=0 moments are exactly one in the saved outputs.

All cutoffs share one table domain, all original mass panels, supplied
power hashes and angular averages. GSL256->512 changes moments/terms
<=3.16e-7 fractionally; internal boost4->8 <=9.46e-6. Their common
errors cancel in cutoff differences: the largest contrast shift is
1.81e-9 fractionally. Original-versus-expanded table domain at the same
1e6 cutoff changes ingredients <=1.23e-5 fractionally; keep this separate.

At masses 1e6,1e4,1e2, respectively, about 5.44%,40.52%,65.07% of sigma^2
comes from above the supplied k endpoint. OneCov's native hmf filter on
the same continued power confirms the 1e5->1e7 integration-tail change
in sigma is <=3.48e-8 fractionally; it agrees with core sigma within
1.93e-6. This validates numerical coverage of the chosen continuation,
not its physical accuracy. Do not describe the lower-cutoff run as a
calibrated improvement to small-halo physics. The subsequent production
adoption is a separate validation, described below.

Halo moments only, 3 redshifts x 41 k x 861 unordered pairs, GL512:
11-repeat mean/stddev at 8 threads are 5.58+/-0.55,7.14+/-0.82,
9.17+/-1.33 ms at the three cutoffs. These exclude first-use tables,
power/tree/response generation and survey projection. No production C/C++
source, installed interface or OneCov numerical code was changed.

Full-matrix follow-up: `run_cutoff_covariance.py` reuses the official LSST
YAML and production backend, including Gaussian non-Limber. Keep AB1,
core boost1 and integration0 (96 nodes), unlike the refined ingredient
test. Complete 1560x1560 G/SSC/cNG/total archives are under
work/cutoff_full_{native6,wide6,wide4,wide2}. The three wide runs share
the sigma domain down to1e2 and identical CAMB inputs. No likelihood mask
is applied. All totals pass Cholesky; minimum normalized eigenvalue is
about2.27e-4. Gaussian components are bitwise equal.

For wide6 versus wide2, largest total variance-mode change is4.5581e-8
fractionally; wide4 versus wide2 gives1.7958e-9. Component differences
against the reference total are also small: SSC2.5426e-8 and cNG4.3651e-8
for6vs2. Changing only the sigma-table domain at Mmin1e6 instead produces
5.9639e-6 in total modes; native6 versus wide2 gives5.9578e-6. Keep that
regridding effect separate from integration-cutoff dependence.
`compare_full_covariance.py` checks matched settings, every matrix entry,
component sums, positivity and generalized modes and draws all four full
difference matrices. Results: results/cutoff_full_*_20261006.json.

Three sequential fresh processes at each cutoff: mean/std total construction
49.57+/-2.12 s,50.88+/-0.83 s,50.18+/-1.69 s for Mmin1e6/1e4/1e2.
Includes first-use tables; excludes CAMB setup and file writing. The
scatter prevents a precise overall slowdown claim; the64% number remains
an ingredient-only result. `report_cutoff_covariance.py` verifies bitwise
repeatability of all physical outputs and inputs, and saves per-stage
and per-run timings. No simultaneous numerical jobs were used.

The1e4 cutoff is a reasonable compromise for reducing completion weight,
not proof of a physically more accurate covariance. The1e6 cutoff already
has negligible measured covariance sensitivity. After reviewing these
results the maintainer authorized adopting1e4 in production and running
all project suites. That separate change must validate its actual sigma
table domain, preserve the original mass-panel nodes, and retain the
additive I11 prescription and multiplicity/bias conventions. The cited
Mead Appendix A prescription is Eq52 specialized to matter, not the
constant additive Eq50; no analogous I12/I13 completion is implemented.

Production adoption check: work/cutoff_full_production4 uses the actual
1e4 table domain and all seven projects' new common mass-panel helper.
The full LSST covariance takes 50.628 s in this single eight-thread run.
It is positive definite, with G bitwise equal to all archived controls.
Maximum total variance-mode changes are 3.13474e-6 against native6,
3.08926e-6 against wide4 and 3.08907e-6 against wide2, fractionally.
All component matrices were compared separately, with no cuts or repairs.
Results and complete difference plots: cutoff_production4_vs_*_20261006.json
and figures/cutoff_production4_vs_*.png. Their source hashes pin the
uncommitted production state as well as the compiled-interface fingerprint.
Repository regression validation is recorded in the core skill's
covariance_mass_cutoff.md. The three-repeat wide-table timing remains the
controlled cutoff-cost comparison; do not turn the single adoption run
into a repeated timing claim.

The diagnostic builder now expects the production 1e4 initializer.
Use --log10-min 6 for an isolated original-domain control. The ingredient
runner restricts its requested mass cutoffs to that build's supported
range, so reproducing the old control no longer relies on an old installed
production library. OneCov itself remains unchanged.
The rebuilt 1e6 control was checked after adoption: its ingredient run
completes, and all four full covariance components reproduce the original
archive bitwise. Local outputs are mass_cutoff_original_reproduced_check
and cutoff_full_native6_reproduced under work/.

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

## Controlled stage timing, 2026-10-06

`time_nongaussian.py` measures resident-input SSC/cNG projection and native
halo-moment calls, without the diagnostic profiler. OneCov projection
statements are extracted by AST from cov_ell_space.py without arithmetic
changes; only the destination assignment becomes a return. Source and
expression hashes are recorded. This is an explicit diagnostic adapter,
not a public full-method timing or a replacement numerical implementation.
Each output reproduces the earlier actual native method export.

Two sequential processes per backend on quiet M2 Pro, OMP8/BLAS1. SSC
100x100, radial300/2401: Cocoa 0.236/1.098 ms, OneCov 9.896/110.437 ms.
cNG8x8, radial300/601: Cocoa 0.210/0.348 ms, OneCov 0.091/0.107 ms.
No mask, response, halo or angular-table generation enters those timers.
Input layout preparation is excluded on both sides; fresh output and
intermediate allocations are included. Report all results, including
OneCov's advantage on the small cNG blocks. Never infer full-survey ratios.

Native I11/I12/both I13/I04 at three z and 200 k nodes: Cocoa 19.93 ms,
OneCov 2497.84 ms. The native fits differ. Cocoa's extra I02 and OneCov's
native ancillary work stay timed. First calls/setup recorded separately;
the headline excludes first-use sigma tables and cosmology. Saved native
I11/diagonal I12 reproduced exactly. Refining mass400->800 vs GL96->256
changes all moments at k=.001..10 by <=1.45e-5 vs1.30e-6; high-k tails
of the requested 1e-5..92.3 grid reach 1.01% vs0.00514%. No blanket
convergence claim over those tails. See nongaussian_timing_20261006.json.

README comparisons of what each code does use explicitly named bullets.
Explain f(nu)dnu as mass fraction, distinguish normalization from quadrature
accuracy, and never interpret Cocoa's imposed int b f=1 as superior
integration. Both native I11 paths add unresolved-low-mass completion;
OneCov alone retains its finite-range bias divisor in I12/I13.

## Sigma and abundance diagnosis, 2026-10-06

`diagnose_sigma.py` exports the actual Cocoa power reader, sigma and mass
slope, holding its input power fixed while refining internal tables 1/2/4.
The OneCov environment calls its unchanged hmf TopHat on those arrays.
This is a two-code shared-input comparison, not a third reference library.
Native hmf sigma is reproduced to roundoff before changing inputs.

At z=.1,.5,1 over 1e10..1e15 Msun/h, native maximum sigma offsets are
0.03238%,0.03044%,0.02959%. Replacing P on the same 200 k nodes changes
sigma by <=0.01608%; dense integration of that same Cocoa reader changes
it by <=0.01594%. Range continuation contributes <=0.000158%, the density
convention <=0.001607%. Same-power/same-radius FFTLog plus table residuals
are <=0.000334%,0.000266%,0.001175%; internal boost4 reduces the largest
to0.000611%. Maxima occur at different masses and must not be added.

OneCov native power200->400->800->1600->3200->6400->12800 was tested;
the last sigma change is1.40e-7 fractionally. Native normalization changes
with this grid too. The initial prepare_lsst_y1.py sigma8 conversion is
a separate CAMB call, not the exact forecast call: supplied0.826717834
vs actual forecast CAMB0.826704159. Forecast nnu3.046 vs pilot3.044.
OneCov's dense filter on forecast CAMB interpolated P gives0.826616612;
on the core reader gives0.826613956. Core FFTLog at R8, internal boost4,
gives0.826615311. Do not blame the resulting amplitude mismatch on FFTLog
or silently overwrite the frozen original inputs. An explicitly matched
future native comparison should resolve this amplitude/radiation setup.
Matching sigma8 for the dense spectra at each z leaves <=0.0081% shape
differences. The detailed source of that remainder is not yet isolated.

Abundance ratio factorization into density, f at matched nu, peak shift
and mass slope reproduces the archived ratio within6.7e-16. Matched-nu
multiplicity ratios are constant over masses:1.00648649,1.02622216,
1.04832034. This explains the growing z1 abundance offset. After dividing
out this amplitude, residual ranges are[-.1908,.1246]%,[-.2163,.1133]%,
[-.2714,.0994]%. It is not a 0.03% sigma error amplified into5%.
`report_sigma.py` writes the quantitative record and science figure;
no numerical source of either code was changed.

Response citations must distinguish implementation from context. Cocoa's
fractional halo response uses corrected Takada-Hu Eq44, then transfers it
to target nonlinear P. OneCov uses the linear slope (no I11 derivative).
Barreira-Krause-Schmidt2018 is broader density/tidal response context;
neither tested isotropic halo prescription implements its full treatment.
