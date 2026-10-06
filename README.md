# CoCoA vs OneCovariance: covariance comparison

This repository will compare covariance calculations from
[CoCoA/CosmoLike](https://github.com/CosmoLike/cocoa) and
[OneCovariance](https://github.com/rreischke/OneCovariance).
It follows the reproducible-study format of
[CCL-benchmark](https://github.com/vivianmiranda/CCL-benchmark): matched
inputs, component comparisons, convergence checks, figures and timing logs.

The aim is to understand where the codes agree, explain differences in
their physical models and numerical methods, and measure execution time
at settings whose numerical accuracy has been checked.

**Status:** study setup. Comparison scripts and measured results will be
added as each stage is completed.

## Contents

1. [Scope](#scope)
2. [Comparison stages](#stages)
3. [Accuracy and execution time](#validation)
4. [Reproducing the comparison](#reproduction)

## Scope <a name="scope"></a>

The initial survey configurations are LSST Y1 and Roman real. We will start
with small parts of their cosmic-shear, galaxy–galaxy-lensing and clustering
covariances. Full survey matrices are a later goal, subject to measured
runtime and memory requirements. Fourier-space comparisons will precede
real-space ones.

The covariance contributions will be kept separate:

| Contribution | Physical origin |
| --- | --- |
| Gaussian (G) | Products of two-point correlations, including shot and shape noise. |
| Super-sample (SSC) | Modes larger than the survey change the fluctuations measured inside it. |
| Connected non-Gaussian (cNG) | Connected four-point correlations beyond the Gaussian contribution. |

The first comparisons will use massless neutrinos, Limber spectra, zero
intrinsic alignment and a simple survey footprint. Galaxy samples will
use matched linear bias rather than OneCovariance's HOD model. Non-Limber
and intrinsic-alignment comparisons will follow separately where the
implemented models can be matched.

The DES cluster 6×2pt + counts calculation is outside the initial scope.
OneCovariance's photometric/spectroscopic 6×2pt example describes a
different set of observables.

## Comparison stages <a name="stages"></a>

1. **Gaussian covariance from shared angular spectra.** Supply identical
   spectra, noise levels and bin definitions. Start with one source bin
   and a few multipole bands to isolate covariance assembly from spectrum
   generation and establish the first runtime measurement.
2. **Matter trispectrum and SSC ingredients.** Compare the 1h, 2h, 3h
   and 4h contributions, including both two-halo partitions, and the
   matter-power response entering SSC. Begin at one redshift with a few
   wavenumbers, including unequal wavenumber pairs. Match halo definitions,
   mass functions, concentration relations and integration ranges first.
3. **Selected survey blocks.** Compute SSC and cNG separately for a small
   shear block, then add representative galaxy and cross-bin blocks.
   Compare their sum with G using the same cosmology, redshift
   distributions, galaxy biases, densities, shape noise and binning.
4. **Diagnose remaining differences.** Vary one choice at a time and
   inspect intermediate quantities. In particular, distinguish Cocoa's
   full-sky real-space transforms from OneCovariance's flat-sky transforms.
5. **Expand only after the pilots.** Add redshifts, wavenumbers and output
   blocks in small increments. Use their measured costs to decide which
   comparisons fit the laptop and which full matrices belong on a server.

Shared-input tests will isolate individual calculations. Separate runs
with each code generating its own ingredients will test the full pipeline.
Diagnostic adaptations will be documented separately from comparisons
using the released implementations.

### Keeping laptop runs small

OneCovariance's runtime on this machine has not yet been measured. Each
pilot will have a stated wall-time budget and record peak memory as well
as elapsed time. Full matrices and broad parameter sweeps will not be
the default laptop workload.

The timing breakdown must separate initialization, shared halo tables and
projection/assembly. A single small block cannot predict the full cost:
some work is paid once per cosmology, while other work grows with the
number of wavenumber pairs, redshift samples or covariance blocks.

Small tests must still sample low and high wavenumbers, more than one
redshift and off-diagonal entries before a broader agreement claim.
For real-space tests, reduce the number of requested angular bins while
keeping the multipole range and integration accuracy needed to converge
those bins. A smaller output matrix is not a reason to lower the physical
integration cutoff.

Reuse saved input tables for projection diagnostics, recording their
provenance. Report these timings separately from runs that generate the
tables themselves.

## Accuracy and execution time <a name="validation"></a>

Each code will first be compared with its own higher-accuracy calculation.
Neither code is assumed to be the reference truth merely because its
settings are more expensive.

Diagnostics will include diagonal variance ratios, correlation matrices,
component differences and positive definiteness of the total covariance
after identical scale cuts. Generalized eigenvalues will quantify changes
in variance across all linear combinations of the measurements. Fisher
comparisons can then assess the impact on constraints for an explicitly
chosen parameter set.

The data-vector criterion $`\Delta\chi^2 < 0.2`$ is not automatically a
covariance acceptance criterion. Numerical convergence and differences
between physical models will be reported separately.

Timings will use Cocoa's production CLI and OneCovariance's normal runner,
with one numerical job at a time on a quiet machine. Reports will record
hardware, thread counts, code revisions and numerical settings. Full
covariance construction will include first-use tables; initialization,
file writing and plotting will be identified separately. Repeated runs
will provide a mean and a measure of timing variation.

## Reproducing the comparison <a name="reproduction"></a>

As the study develops, this repository will contain:

| Planned directory | Contents |
| --- | --- |
| `scripts/` | Python runners, comparison diagnostics and plotting scripts. |
| `configs/` | Survey settings and numerical configurations for both codes. |
| `inputs/` | Shared input tables and records of how they were generated. |
| `results/` | Comparison outputs, machine-readable metrics and timing logs. |
| `figures/` | Figures generated from the saved results. |

Every reported result will identify the code revisions, configuration,
input provenance and command needed to reproduce it. The README will
explain the physics and summarize the findings, with links to the saved
evidence and any unresolved differences.
