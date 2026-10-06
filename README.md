# CoCoA vs OneCovariance: covariance comparison

This repository compares covariance calculations from
[CoCoA/CosmoLike](https://github.com/CosmoLike/cocoa) and
[OneCovariance](https://github.com/rreischke/OneCovariance).
It follows the reproducible-study format of
[CCL-benchmark](https://github.com/vivianmiranda/CCL-benchmark): matched
inputs, component comparisons, convergence checks, figures and timing logs.

The aim is to understand where the codes agree, explain differences in
their physical models and numerical methods, and measure execution time
at settings whose numerical accuracy has been checked.

This is an **accuracy comparison first**. We are not optimizing or
rewriting OneCovariance. Timing comparisons follow once the physical
choices and numerical accuracy of both calculations are understood.

**First comparison completed:** Gaussian covariance assembly agrees for
shared LSST Y1 angular spectra, matched noise and matched band weights.
Three small tests cover one source bin alone and two lens bins with that
source, including all crossed spectra. Both codes give positive definite
total matrices in every case.

| Shared-spectrum test | Matrix | Integer multipoles | Largest variance-mode difference |
| --- | --- | --- | ---: |
| Shear, source bin 3 | 5×5 | 30–149 | 0.000017% |
| 3×2pt, lens bins 1 and 2, source bin 3 | 30×30 | 30–149 | 0.000465% |
| Same 3×2pt subset | 30×30 | 1500–1619 | 0.000046% |

The differences are consistent with **OneCov's text-output rounding**:
OneCov saves totals at seven significant digits and split components at
five significant digits. The mode comparison considers every linear
combination of the selected bandpowers, not just diagonal variances:
it reports the largest $`|\lambda-1|`$ for
$`C_{\mathrm{OneCov}}v=\lambda C_{\mathrm{Cocoa}}v`$, expressed as a percentage.

See the [comparison record](results/gaussian_assembly_20261005.json) and
[reproduction steps](#matched-gaussian). This establishes the Gaussian
contractions, noise normalization and binning for **shared spectra**.
It does not establish agreement of native spectra or halo ingredients.
The [SSC comparison](#ssc-comparison) and
[cNG projection comparison](#connected-comparison) separately test the
non-Gaussian contributions with shared inputs. The
[halo study](#halo-comparison) identifies native-model differences, and the
[trispectrum study](#trispectrum-comparison) isolates an off-diagonal
2-halo assembly discrepancy. The timing comparison below covers Gaussian
components and band averaging only.

## Contents

1. [Scope](#scope)
2. [Comparison stages](#stages)
3. [Accuracy and execution time](#validation)
4. [Installation and compilation](#installation)
5. [Reproducing the comparison](#reproduction)
6. [Matched Gaussian assembly](#matched-gaussian)
7. [SSC comparison](#ssc-comparison)
8. [Halo-model ingredients](#halo-comparison)
9. [Separated halo trispectra](#trispectrum-comparison)
10. [Connected non-Gaussian projection](#connected-comparison)

## Scope <a name="scope"></a>

**LSST Y1 is the benchmark.** We start with small parts of its cosmic-shear,
galaxy–galaxy-lensing and clustering covariances. Full survey matrices are
a later goal, subject to measured runtime and memory requirements.
Fourier-space comparisons precede real-space ones. Other surveys can be
added after this comparison is understood.

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

## Gaussian results and execution time <a name="validation"></a>

The correlation matrices agree at both multipole ranges. The right panels
magnify the saved-matrix differences by a million; their small residuals
come from OneCov's seven-significant-digit text output.

![Gaussian correlations and residuals](results/figures/gaussian_matrices.png)

The variance components also agree. Sample variance dominates the low-ell
galaxy spectra; shape noise dominates the high-ell shear variance for these
LSST Y1 densities. Lines show Cocoa and open circles show OneCov. Each
group contains five bands, with multipole increasing from left to right.

![Gaussian variance components](results/figures/gaussian_components.png)

Vector figures: [matrices](results/figures/gaussian_matrices.pdf),
[components](results/figures/gaussian_components.pdf).

### Gaussian timing differences

**Apple M2 Pro, macOS 13.7.5, eight OpenMP threads.** Both codes receive
the same spectra already in memory and compute all three Gaussian parts.
Cocoa uses its production `_interface`; OneCov uses its unmodified
`covELL_gaussian` method.

| Gaussian case | CoCoA (ms) | OneCov (ms) | OneCov / CoCoA |
| --- | ---: | ---: | ---: |
| Shear, 30 ≤ ℓ < 150 | 0.374 ± 0.091 | 0.786 ± 0.021 | 2.1 |
| 3×2pt, 30 ≤ ℓ < 150 | 0.242 ± 0.014 | 9.647 ± 0.193 | 39.8 |
| 3×2pt, 1500 ≤ ℓ < 1620 | 0.243 ± 0.003 | 9.894 ± 0.202 | 40.7 |

These are repeated-batch means; ± shows the scatter between batches.
The tiny shear case is more variable: Cocoa's two run means were 0.327 and
0.422 ms. The [timing record](results/gaussian_timing_20261005.json) retains
all samples, first-call times, setup times and numerical checks.

The timer includes numerical output allocation, but excludes spectrum
generation, initialization and file writing. OneCov's final rearrangement
of its blocks into one matrix is also outside the timer. These are
**small Gaussian component timings, not full-survey covariance runtimes**.

## Installation and compilation <a name="installation"></a>

We use Cocoa's Python 3.11 Conda base and a private `.local` environment
inside this repository. Follow the
[Cocoa Conda installation recipe](https://github.com/CosmoLike/cocoa#required_packages_conda)
if that base is not installed. Cocoa and OneCovariance are separate sibling
checkouts; Cocoa must already have CAMB compiled.

Open a fresh Bash terminal in `OneCov-benchmark-/`. Keep the Cocoa export
session in a separate terminal.

**Step :one:**: activate the Conda base.

```bash
conda activate cocoa
```

**Step :two:**: inspect [set_installation_options.sh](set_installation_options.sh).
It specifies the code paths, studied revisions and Python package versions.

**Step :three:**: prepare the private environment and download dependencies.

```bash
source setup_onecov.sh
```

**Step :four:**: compile the native extensions from local sources.

```bash
source compile_onecov.sh
```

**Step :five:**: activate the comparison environment.

```bash
source start_onecov.sh
```

| File | Purpose |
| --- | --- |
| [set_installation_options.sh](set_installation_options.sh) | Select code paths and pinned package versions. |
| [setup_onecov.sh](setup_onecov.sh) | Create `.local`, install Python dependencies and download healpy sources. |
| [compile_onecov.sh](compile_onecov.sh) | Build healpy and OneCovariance's bundled Levin extension offline; check imports. |
| [start_onecov.sh](start_onecov.sh) | Activate the private environment. |
| [stop_onecov.sh](stop_onecov.sh) | Restore the shell's previous Python environment. |

The scripts reuse the existing OneCovariance checkout and Cocoa's CAMB.
They do not modify either numerical code or install packages into Cocoa.
The [environment notes](docs/environment.md) describe the shared libraries
and the initial pilot environment.

### Starting and stopping later sessions

After installation, open a fresh Bash terminal in `OneCov-benchmark-/`.
Setup and compilation do not need to be repeated for each calculation.

**Step :one:**: activate the Conda base.

```bash
conda activate cocoa
```

**Step :two:**: activate the comparison's private `.local` environment.

```bash
source start_onecov.sh
```

**Step :three:**: after the calculation, leave the private environment.

```bash
source stop_onecov.sh
```

## Reproducing the comparison <a name="reproduction"></a>

### What the pilot represents

The exporter reads Cocoa's `projects/lsst_y1/covariance` configuration.
The default subset is source bin 3 and lens bins 1 and 2, retaining their
original LSST Y1 identities and densities.

| Input | LSST Y1 forecast choice |
| --- | --- |
| Area | 12,300 deg² |
| Source density | 2 arcmin⁻² in the selected source bin |
| Lens density | 3.6 arcmin⁻² in each selected lens bin |
| Shape dispersion | 0.26 per component |
| Lens biases | 1.72716 and 1.65168 for lens bins 1 and 2 |
| Cosmology | Ωm = 0.3, Ωb = 0.05, h = 0.7, ns = 0.965, As = 2.1×10⁻⁹ |
| First output | Five Fourier bands between ℓ = 30 and 3,000 |

Selecting one source bin does **not** put all five bins' galaxies into it.
The scripts preserve its density and export its actual n(z) column. These
are the project's forecast choices, not its supplied likelihood covariance.

### Exporting the Cocoa inputs

Follow the [main Cocoa setup](https://github.com/CosmoLike/cocoa#cobaya_base_code_examples)
and [LSST Y1 covariance instructions](https://github.com/CosmoLike/cocoa_lsst_y1#computing_covariances).
We assume Cocoa and LSST Y1 are installed, the shell is Bash, and
`cocoa/`, `OneCovariance/` and `OneCov-benchmark-/` are sibling directories.

Open a second, fresh Bash terminal in `OneCov-benchmark-/`.

**Step :one:**: activate the Conda base.

```bash
conda activate cocoa
```

**Step :two:**: enter Cocoa's runtime directory.

```bash
cd ../cocoa/Cocoa
```

**Step :three:**: activate Cocoa's private environment.

```bash
source start_cocoa.sh
```

**Step :four:**: enable the LSST Y1 project and covariance bindings.

```bash
unset IGNORE_COSMOLIKE_LSST_Y1_CODE IGNORE_COSMOLIKE_LSST_Y1_COVARIANCE
```

**Step :five:**: compile the project interface.

```bash
source ./projects/lsst_y1/scripts/compile_lsst_y1.sh
```

**Step :six:**: select eight OpenMP threads using your platform's settings.

- Linux

  ```bash
  export OMP_NUM_THREADS=8 OMP_PROC_BIND=close \
    OMP_PLACES=cores OMP_DYNAMIC=FALSE \
    OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ```

- macOS (arm)

  ```bash
  export OMP_NUM_THREADS=8 OMP_PROC_BIND=disabled \
    OMP_PLACES=cores OMP_DYNAMIC=FALSE \
    OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ```

**Step :seven:**: return to the benchmark directory.

```bash
cd ../../OneCov-benchmark-
```

**Step :eight:**: export the LSST Y1 inputs.

```bash
python scripts/prepare_lsst_y1.py \
  --cocoa ../cocoa/Cocoa --output work/lsst_y1
```

This writes the selected distributions, constant bias tables, shared
noise-free angular spectra, CAMB nonlinear power and a manifest containing
revisions and input hashes. It computes sigma8 from the specified As for
OneCovariance's amplitude input. It does not compute a covariance matrix.

Use `--source-bin 4 --lens-bins 2 3` to select other original LSST bins.
The two lens bins must be distinct. At the studied OneCov revision, its
supplied-bias path fails for a single lens bin because the reader and
constructor disagree about the array shape. Two real lens bins avoid that
problem without modifying OneCovariance.

### Running the small Gaussian case

We assume installation is complete and the inputs are in `work/lsst_y1`.
Open a fresh Bash terminal in `OneCov-benchmark-/`. Run one numerical job
at a time.

**Step :one:**: activate the Conda base.

```bash
conda activate cocoa
```

**Step :two:**: activate the comparison's private environment.

```bash
source start_onecov.sh
```

**Step :three:**: select eight OpenMP threads using your platform's settings.

- Linux

  ```bash
  export OMP_NUM_THREADS=8 OMP_PROC_BIND=close \
    OMP_PLACES=cores OMP_DYNAMIC=FALSE \
    OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ```

- macOS (arm)

  ```bash
  export OMP_NUM_THREADS=8 OMP_PROC_BIND=disabled \
    OMP_PLACES=cores OMP_DYNAMIC=FALSE \
    OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ```

**Step :four:**: compute the five-band shear Gaussian covariance.

```bash
python scripts/run_onecov.py \
  --onecov ../OneCovariance --inputs work/lsst_y1 \
  --output work/shear_gaussian --timeout 180
```

**Step :five:**: check the saved shear covariance.

```bash
python scripts/check_result.py work/shear_gaussian
```

The default uses shared Cocoa angular spectra and produces a **5×5 shear
Gaussian covariance**. `check_result.py` verifies finite entries, symmetry,
and positive definiteness.

**Step :six:**: include galaxy clustering, galaxy–shear and their cross blocks.

```bash
python scripts/run_onecov.py \
  --inputs work/lsst_y1 --case 3x2 \
  --output work/small_3x2_gaussian --timeout 180
```

**Step :seven:**: check the saved 3×2pt covariance.

```bash
python scripts/check_result.py work/small_3x2_gaussian
```

This gives a **30×30 matrix**, including the cross spectrum between the two
lens bins. The check covers matrix finiteness, symmetry and positivity.

> [!NOTE]
> OneCovariance averages integer multipoles uniformly in these Gaussian
> bands. Cocoa's Fourier operator uses mode-count weights. Identical band
> edges alone therefore do not define identical estimators. The comparison
> below supplies uniform weights to Cocoa's existing production Gaussian
> kernel. This benchmark adapter does not change Cocoa's default estimator.

### Comparing matched Gaussian assembly <a name="matched-gaussian"></a>

Use the two terminals prepared above: one with Cocoa active, the other
with `start_onecov.sh` active. Both should be in `OneCov-benchmark-/`, with
the platform's eight-thread settings. Run the following steps sequentially.

These tests use five bands in a narrow multipole interval. Exporting each
integer multipole removes interpolation from the comparison. This matters:
on a coarse grid, interpolating a product of spectra is different from
interpolating each spectrum and then multiplying.

**Step :one:**: in the **Cocoa terminal**, export the low-multipole inputs.

```bash
python scripts/prepare_lsst_y1.py \
  --cocoa ../cocoa/Cocoa --output work/lsst_y1_integer_30_150 \
  --integer-ell-range 30 150
```

The export includes both endpoints. The covariance sums include the lower
band edge and exclude the upper edge, so these bands use multipoles 30–149.

**Step :two:**: in the **OneCov terminal**, compute the shear covariance.

```bash
python scripts/run_onecov.py \
  --inputs work/lsst_y1_integer_30_150 --case shear \
  --band-limits 30 150 --output work/assembly_shear_30_150
```

**Step :three:**: in the **OneCov terminal**, compute all small 3×2pt blocks.

```bash
python scripts/run_onecov.py \
  --inputs work/lsst_y1_integer_30_150 --case 3x2 \
  --band-limits 30 150 --output work/assembly_3x2_30_150
```

**Step :four:**: in the **Cocoa terminal**, compare the shear calculation.

```bash
python scripts/compare_gaussian.py work/assembly_shear_30_150 \
  --cocoa ../cocoa/Cocoa --output work/reviewed_gaussian/shear_30_150
```

**Step :five:**: in the **Cocoa terminal**, compare the 3×2pt calculation.

```bash
python scripts/compare_gaussian.py work/assembly_3x2_30_150 \
  --cocoa ../cocoa/Cocoa --output work/reviewed_gaussian/3x2_30_150
```

**Step :six:**: in the **Cocoa terminal**, export the higher-multipole inputs.

```bash
python scripts/prepare_lsst_y1.py \
  --cocoa ../cocoa/Cocoa --output work/lsst_y1_integer_1500_1620 \
  --integer-ell-range 1500 1620
```

**Step :seven:**: in the **OneCov terminal**, compute this 3×2pt covariance.

```bash
python scripts/run_onecov.py \
  --inputs work/lsst_y1_integer_1500_1620 --case 3x2 \
  --band-limits 1500 1620 --output work/assembly_3x2_1500_1620
```

**Step :eight:**: in the **Cocoa terminal**, compare the higher multipoles.

```bash
python scripts/compare_gaussian.py work/assembly_3x2_1500_1620 \
  --cocoa ../cocoa/Cocoa --output work/reviewed_gaussian/3x2_1500_1620
```

`compare_gaussian.py` calls Cocoa's production `_interface` with the shared
spectra, diagonal shot/shape noise and uniform band weights, then compares
the resulting components with OneCov's output.
No numerical source in either code is modified.

Each comparison saves `comparison.json` and `matrices.npz`. These retain
sample variance, mixed signal/noise, pure noise and their total separately.
The script checks positive total matrices and generalized variance ratios
between the codes. It fails if agreement exceeds the native text precision.

For the 3×2pt tests, matrix order is $`g_1g_1`$, $`g_1g_2`$, $`g_2g_2`$,
$`g_1\gamma`$, $`g_2\gamma`$, $`\gamma\gamma`$, with five bands per spectrum.
Including $`g_1g_2`$ checks cross-bin contractions even though this spectrum
need not belong to the likelihood's data vector.

Agreement of angular spectra generated separately by each code remains
a separate test. Shared CAMB power can isolate their projection before
native power and its numerical convergence are compared.

### Reproducing Gaussian plots and timings

Use the prepared Cocoa and OneCov terminals from the preceding section.
The commands below time the low-ell 3×2pt case. To time the other cases,
replace `3x2_30_150` with `shear_30_150` or `3x2_1500_1620` in both paths.

**Step :one:**: in the **OneCov terminal**, time its Gaussian components.

```bash
python scripts/time_gaussian.py work/reviewed_gaussian/3x2_30_150 \
  --backend onecov --output work/gaussian_timing/onecov_3x2_30_150
```

**Step :two:**: after it finishes, in the **Cocoa terminal**, time Cocoa.

```bash
python scripts/time_gaussian.py work/reviewed_gaussian/3x2_30_150 \
  --backend cocoa --output work/gaussian_timing/cocoa_3x2_30_150
```

Each run saves `timing.json`, including 31 batches and a separately
measured first call. Cocoa batches contain 100 calls; OneCov batches
contain 10. Both recompute outputs on every call. Repeats require a new
output directory, preserving earlier measurements.

**Step :three:**: in either terminal, regenerate the figures from the
validated matrices.

```bash
python scripts/plot_gaussian.py --comparisons work/reviewed_gaussian \
  --output results/figures
```

This writes PNG and PDF versions of the matrix and component figures.
Execution times are reported in the table above.

## SSC comparison <a name="ssc-comparison"></a>

The projected comparison covers a **100×100 shear SSC matrix** for LSST
Y1 source bin 3, at multipoles from 30 to 3000. CoCoA and OneCov use the
same matter responses, linear power, radial samples and lensing window.
The matrix includes off-diagonal multipole correlations. This isolates
projection from differences in the two codes' halo prescriptions.

The codes agree to floating-point precision. Replacing only OneCov's
pixelized spherical-cap footprint with CoCoA's analytic cap changes the
matrix by at most **0.040%**, at the same area and mask multipole cutoff.
The metric is $`|\Delta C_{ij}|/\sqrt{C_{ii}C_{jj}}`$; no diagonal-only
agreement is assumed. See the
[SSC projection record](results/ssc_projection_20261005.json).

![SSC projection comparison](results/figures/ssc_projection.png)

[Vector figure](results/figures/ssc_projection.pdf).

Refining the shared calculation gives the following largest matrix
changes, using the finer matrix's diagonal for normalization:

| Refinement | Largest SSC change |
| --- | ---: |
| Radial nodes: 300 → 601 | 0.0430% |
| Radial nodes: 601 → 1201 | 0.0412% |
| Response redshift spacing: 0.1 → 0.05, at 1201 radial nodes | 0.1831% |
| Radial nodes: 1201 → 2401, at spacing 0.05 | 0.0040% |
| Response redshift spacing: 0.05 → 0.025, at 2401 radial nodes | 0.0427% |

Every shared-input CoCoA–OneCov comparison agrees within
$`7\times10^{-15}`$ in the same metric, and all six SSC matrices are
positive definite. These tests establish the shared-input projection;
halo mass/power sampling and native response choices still require their
own comparisons. Small SSC entries and Fisher information are not
certified by the diagonal-normalized convergence metric alone.

### Matter-power response

The first SSC ingredient comparison now checks the **matter-power
response**, $`D=\partial P/\partial\delta_b`$: how power changes inside a
large-scale background overdensity. With shared halo ingredients and a
matched response prescription, Cocoa and OneCov agree to **4×10⁻¹⁶** at
redshifts 0, 0.5 and 1. See the
[response comparison record](results/ssc_response_20261005.json).

This checks the response formula, not the projected SSC covariance.
The default prescriptions differ: OneCov uses the linear-power slope;
Cocoa uses the two-halo slope and transfers the fractional halo response
to its nonlinear power. The two-halo prescription follows the corrected
Eq. 44 of [Takada & Hu (2013)](https://arxiv.org/html/1302.6994v3).
The broader response treatment is described by
[Barreira, Krause & Schmidt (2018)](https://arxiv.org/abs/1711.07467).

On shared ingredients, changing only that slope changes the sampled
responses by at most 0.23%, 0.18% and 0.14%, respectively, over
$`0.001\leq k\leq10\,h/\mathrm{Mpc}`$. This diagnostic excludes nonlinear
power transfer and does not establish convergence of the native forecasts.

To reproduce this ingredient check after the Gaussian example:

**Step :one:**: in the **OneCov terminal**, export the native responses.

```bash
python scripts/ssc_response.py export work/assembly_shear_30_150/onecov.ini \
  --output work/ssc_response_onecov
```

**Step :two:**: in the **Cocoa terminal**, compare the response formula.

```bash
python scripts/ssc_response.py compare work/ssc_response_onecov \
  --output work/ssc_response_comparison
```

## Halo-model ingredients <a name="halo-comparison"></a>

The native halo prescriptions **do not produce identical inputs**. We
compare $`z=0.1,0.5,1`$, masses $`10^{10}`$–$`10^{15}\,M_\odot/h`$, and
wavenumbers $`0.001`$–$`10\,h/\mathrm{Mpc}`$. The following are the largest
absolute differences in CoCoA/OneCov minus one over those sampled ranges:

| Quantity | z = 0.1 | z = 0.5 | z = 1 |
| --- | ---: | ---: | ---: |
| Mass rms fluctuation, sigma(M) | 0.033% | 0.031% | 0.030% |
| Halo abundance, dn/dlnM | 0.8% | 2.8% | 5.0% |
| Native halo bias | 22.8% | 28.9% | 35.9% |
| Native matter-power response | 14.9% | 21.2% | 25.9% |

The raw [Tinker et al. (2010)](https://arxiv.org/abs/1001.3162) bias
formulas agree within $`6\times10^{-16}`$ when evaluated at the same
peak height. OneCov then divides the bias by a finite-mass-range
normalization, measured here as 0.7725, 0.7119 and 0.6415. CoCoA instead
sets the multiplicity normalization through its bias-consistency integral.
These choices explain the large bias offset; it is not a disagreement
in the underlying bias formula.

The concentration choices also differ: OneCov uses
[Duffy et al. (2008)](https://arxiv.org/abs/0804.2486), while CoCoA uses
[Bhattacharya et al. (2013)](https://arxiv.org/abs/1112.5479).
At the same supplied concentration, the NFW profiles differ by less than
$`1.4\times10^{-5}`$ absolutely over the exported grid. The native moment
mass limits and SSC response transfer remain different, as described above.

![Native halo ingredients](results/figures/halo_ingredients.png)

[Vector figure](results/figures/halo_ingredients.pdf) ·
[Detailed results](results/halo_ingredients_20261005.json).

Mass-integration refinement is much smaller than these native differences.
For I11, I02, I12 and the response, OneCov's 400 → 800 mass-node change is
at most **0.0012%**; CoCoA's 96 → 256 nodes per mass panel changes them by
at most **0.000040%**. This checks the mass integration at the supplied
power resolution; it does not certify either prescription's physical accuracy.

**Step :one:**: in the **OneCov terminal**, export the halo quantities.

```bash
python scripts/compare_halo.py export work/shear_ssc/onecov.ini \
  --mass-nodes 800 --output work/halo_onecov_800
```

**Step :two:**: in the **Cocoa terminal**, evaluate the corresponding
CoCoA quantities at the same redshifts, masses and wavenumbers.

```bash
python scripts/compare_halo.py compare work/halo_onecov_800 \
  --integration-accuracy 2 --output work/halo_cocoa_matched_800
```

**Step :three:**: in either terminal, make the halo comparison figure.

```bash
python scripts/plot_halo.py work/halo_onecov_800 work/halo_cocoa_matched_800 \
  --output results/figures
```

### What the bias normalization changes

The Tinker fit is a calibrated halo bias. Its consistency condition refers
to the **full** mass distribution, not an arbitrary numerical mass range:
$`\int b(\nu)f(\nu)\,d\nu=1`$.
[Tinker et al. (2010), Eq. 7](https://arxiv.org/html/1001.3162)

Evaluating the two codes' fitted functions over an extended peak-height
range gives:

| Redshift | OneCov raw $`\int bf\,d\nu`$, extended range | OneCov finite-range divisor | Bias multiplier $`1/N`$ | CoCoA $`\int bf\,d\nu`$ |
| ---: | ---: | ---: | ---: | ---: |
| 0.1 | 0.993555 | 0.772506 | 1.29449 | 1.000000 |
| 0.5 | 0.974448 | 0.711928 | 1.40464 | 1.000000 |
| 1.0 | 0.953907 | 0.641533 | 1.55877 | 1.000000 |

Most of the finite-range deficit therefore comes from excluded low-peak
halos. Dividing the fitted bias by that deficit raises it at **every**
resolved mass. This is an extra prescription; it is not required by the
calibrated Tinker relation.

In the tested OneCov revision, I11 cancels this bias divisor and adds the
missing low-mass contribution explicitly. I12 and I13 retain the divisor.
This different treatment matters for SSC and the multi-halo trispectrum.
At fixed other ingredients, it multiplies 2h(1+3) and 3h by $`1/N`$,
and 2h(2+2) by $`1/N^2`$. At $`z=1`$ these factors are 1.56 and 2.43;
1h and 4h do not receive this correction.

CoCoA instead retains the fitted bias and adjusts the multiplicity
amplitude. Its mass-only integral $`\int f\,d\nu`$ is then 1.00649,
1.02622 and 1.04832; it does **not** simultaneously enforce exact mass
normalization. These are different modeling choices, not interchangeable
implementations of one fit. Integrating extrapolated fits is a consistency
diagnostic, not evidence that the fits are calibrated at arbitrarily low mass.

The [normalization record](results/bias_normalization_20261006.json)
includes the doubled-grid check. To reproduce it:

**Step :one:**: in the **OneCov terminal**, evaluate its native fits.

```bash
python scripts/diagnose_bias.py onecov --config work/shear_ssc/onecov.ini \
  --nodes 8193 --output work/bias_onecov_fine
```

**Step :two:**: in the **Cocoa terminal**, evaluate CoCoA's native fits.

```bash
python scripts/diagnose_bias.py cocoa --nodes 8193 \
  --output work/bias_cocoa_fine
```

## Separated halo trispectra <a name="trispectrum-comparison"></a>

We test 1h, 2h, 3h and 4h at redshifts 0.1, 0.5 and 1, using nine
wavenumbers from 0.001 to 10 h/Mpc. All 45 unordered pairs are included.
The two codes retain their native halo prescriptions for the solid/dashed
curves below; these are **not matched-model predictions**.

![Native halo terms and shared-input assembly](results/figures/trispectrum_terms.png)

With **identical halo moments and angular averages**, the assembly test
finds:

| Contribution | Largest fractional discrepancy |
| --- | ---: |
| 1h | 0 |
| 2h, diagonal pairs | 2.3 × 10⁻¹⁶ |
| 2h, all pairs | **0.1364** |
| 3h | 2.3 × 10⁻¹⁶ |
| 4h | 2.3 × 10⁻¹⁶ |

All entries are fractions; 0.1364 corresponds to **13.64%**. The 1h row
checks that the same supplied one-halo term is copied, not agreement
between independently computed halo integrals.

The off-diagonal 2h discrepancy is localized to the 1+3 halo partitions.
For a covariance configuration $`(K,-K,Q,-Q)`$, their sum is

$$
T^{2h}_{13}=2P_L(K)I^1_1(K)I^1_3(K,Q,Q)
          +2P_L(Q)I^1_1(Q)I^1_3(K,K,Q).
$$

The two moments are different when $`K\ne Q`$.
CoCoA uses both. The sampled OneCov revision uses $`I^1_3(K,Q,Q)`$ in
both terms before mirroring the matrix. Repeating that choice only in
CoCoA's **supplied diagnostic inputs** reduces the discrepancy below
$`3.4\times10^{-16}`$. Neither source implementation was changed.
The partition structure follows [Takada & Hu (2013), Eq. 29](https://arxiv.org/html/1302.6994v3).

Angular averages for unequal wavenumbers agree within
$`2.8\times10^{-7}`$ fractionally when both receive the same linear
power interpolation. Equal pairs approach zero internal wavenumber and
are sensitive to OneCov's corner cutoff and extrapolation. Tightening
its two corner controls from $`10^{-3}`$ to $`10^{-5}`$, then to
$`10^{-7}`$, still changes some native 2h/3h entries substantially.
**The native diagonal angular calculation is not converged by this test.**

By contrast, changing CoCoA's mass/angular rules from 96 to 256 nodes
changes every sampled native term by less than **0.00040%**. Doubling
OneCov's mass grid from 400 to 800 changes them by at most **0.0031%**.
These integration checks do not remove the distinct bias normalization,
concentration relation or OneCov's low-k one-halo damping.

The [trispectrum record](results/trispectrum_20261006.json) separates
shared assembly, native models and refinement diagnostics.

**Step :one:**: in the **OneCov terminal**, export the separated terms.

```bash
python scripts/compare_trispectrum.py export work/shear_ssc/onecov.ini \
  --mass-nodes 800 --output work/trispectrum_800
```

**Step :two:**: in the **Cocoa terminal**, run both comparison scopes.

```bash
python scripts/compare_trispectrum.py compare work/trispectrum_800 \
  --integration-accuracy 2 --output work/trispectrum_final_256
```

## Connected non-Gaussian projection <a name="connected-comparison"></a>

The projection test supplies the same matter trispectrum, lensing window
and radial weights to both codes. It covers eight multipoles between
30 and 3000, including all 64 diagonal and off-diagonal entries.
It tests radial assembly separately from the halo-model differences above.

![Shared connected projection](results/figures/connected_projection.png)

Across ten configurations, the largest
$`|\Delta C_{ij}|/\sqrt{C_{ii}C_{jj}}`$ is
**$`2.3\times10^{-15}`$**. This is agreement of the projection with
shared inputs, not agreement of independently generated halo models.

We also refine OneCov's native **five-band G+cNG** calculation, changing
one control at a time. All ten total matrices are positive definite.
Their correlation matrices have minimum eigenvalues above 0.957.

| Refinement | Largest cNG change, variance-scaled | Largest G+cNG variance-mode change |
| --- | ---: | ---: |
| Wavenumber nodes: 9 → 17 | 38.29% | 2.008% |
| Wavenumber nodes: 17 → 33 | 13.11% | 0.506% |
| Wavenumber nodes: 33 → 65 | 9.59% | 0.276% |
| Wavenumber nodes: 65 → 129 | 2.97% | 0.0903% |
| Radial nodes: 300 → 601 | 0.0375% | 0.0286% |
| Trispectrum redshift step: 0.5 → 0.25 | 7.08% | 0.645% |
| Trispectrum redshift step: 0.25 → 0.125 | 3.15% | 0.311% |
| Mass nodes: 400 → 800 | 0.00130% | 0.000128% |
| Corner controls: $`10^{-3}`$ → $`10^{-5}`$ | 2.95% | 0.0884% |

The cNG column divides each entry change by the finer cNG diagonal rms
product. The total column tests every linear combination of the five
bandpowers using generalized eigenvalues. These two diagnostics need
not be similar: Gaussian noise can make a sizable component change small
in the total covariance.

**The coarse pilot is not a converged native cNG prediction.** In
particular, the final redshift refinement still exceeds a 0.1% total-mode
criterion. These tests do not establish full LSST or Fisher convergence,
nor do they remove the native-model and 2h partition differences above.
The [complete record](results/connected_20261006.json) gives settings,
input hashes and every comparison.

**Step :one:**: in the **OneCov terminal**, compute the native shear block.

```bash
python scripts/compare_connected.py export work/shear_ssc/onecov.ini \
  --k-nodes 65 --radial-nodes 601 --delta-z 0.25 \
  --output work/connected_redshift
```

**Step :two:**: in the **Cocoa terminal**, project its supplied matter tables.

```bash
python scripts/compare_connected.py compare work/connected_redshift \
  --output work/connected_cocoa_redshift
```

**Step :three:**: in either terminal, draw the trispectrum and cNG figures.

```bash
python scripts/plot_connected.py work/trispectrum_final_256 \
  work/connected_cocoa_redshift --output results/figures
```

### Running the projected non-Gaussian pilots

We assume installation is complete and the inputs are in `work/lsst_y1`.
Open a fresh Bash terminal in `OneCov-benchmark-/`. Run the two calculations
sequentially.

**Step :one:**: activate the Conda base.

```bash
conda activate cocoa
```

**Step :two:**: activate the comparison's private environment.

```bash
source start_onecov.sh
```

**Step :three:**: select eight OpenMP threads using your platform's settings.

- Linux

  ```bash
  export OMP_NUM_THREADS=8 OMP_PROC_BIND=close \
    OMP_PLACES=cores OMP_DYNAMIC=FALSE \
    OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ```

- macOS (arm)

  ```bash
  export OMP_NUM_THREADS=8 OMP_PROC_BIND=disabled \
    OMP_PLACES=cores OMP_DYNAMIC=FALSE \
    OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ```

**Step :four:**: compute the shear Gaussian and SSC contributions.

```bash
python scripts/run_onecov.py \
  --inputs work/lsst_y1 --terms ssc --spectra native \
  --output work/shear_ssc --timeout 600
```

**Step :five:**: check the saved Gaussian plus SSC covariance.

```bash
python scripts/check_result.py work/shear_ssc
```

### Reproducing the SSC matrix comparison

After the shear SSC pilot above, use the two prepared environments.

**Step :one:**: in the **OneCov terminal**, export its SSC matrix and
the inputs used by its projection.

```bash
python scripts/compare_ssc.py export work/shear_ssc/onecov.ini \
  --output work/ssc_projection_300
```

**Step :two:**: in the **Cocoa terminal**, compute the same matrix with
CoCoA's production kernels and compare every entry.

```bash
python scripts/compare_ssc.py compare work/ssc_projection_300 \
  --output work/ssc_comparison_300
```

**Step :three:**: in either terminal, make the covariance comparison plot.

```bash
python scripts/plot_ssc.py work/ssc_comparison_300 --output results/figures
```

Use `--radial-nodes`, `--delta-z` and `--mass-nodes` on the export to
refine one numerical input at a time, saving each run in a new directory.
The comparison uses OneCov's actual radial integration rule in both codes.
Native halo inputs and response prescriptions are separate comparisons.

### Running a connected non-Gaussian pilot

**Step :one:**: in the **OneCov terminal**, compute the shear Gaussian and
connected contributions.

```bash
python scripts/run_onecov.py \
  --inputs work/lsst_y1 --terms connected --spectra native \
  --output work/shear_connected --timeout 600
```

**Step :two:**: check the saved Gaussian plus connected covariance.

```bash
python scripts/check_result.py work/shear_connected
```

`--terms ssc` computes G+SSC; `--terms connected` computes G+cNG.
The native list and separate matrix files preserve the contributions.
The pilot's eight-point trispectrum table measures feasibility, **not
converged cNG accuracy**. Halo prescriptions and footprint normalization
still need to be matched before interpreting differences from Cocoa.

The three spectrum choices isolate different calculations:

| `--spectra` | What OneCovariance receives |
| --- | --- |
| `shared-cells` | Cocoa angular spectra; isolates Gaussian covariance assembly. |
| `shared-power` | Cocoa CAMB nonlinear P(k,z); OneCov performs the angular projection. |
| `native` | Cosmology and survey inputs; OneCov generates its own spectra. |

Supplied power does not replace every linear-power calculation inside the
halo response or trispectrum. Shared-input and native runs must retain
their separate labels.

### Inspecting and refining a run

Every run writes `onecov.ini`, `input_manifest.json`, `preflight.log`,
`native.log` and `run.json`. The latter records status, revisions, package
versions, worker count, wall time and peak child memory. A timeout stops
the child process group and keeps partial output. Retry in a new directory.

The recorded wall time covers the fresh native CLI, including imports,
setup, computation and writing. It is a feasibility measurement, not yet
the separated timing breakdown needed for a performance comparison.
Peak child memory is not a sum over simultaneous worker processes.

Copy [onecov_pilot.ini](configs/onecov_pilot.ini), change one numerical
control and pass it with `--template`. `--prepare-only` writes an INI without
importing or running OneCovariance. Refining the internal spectrum grid
does not refine a supplied C_ell table; regenerate that input separately.

Generated files live under ignored `work/`. Small reviewed validation
summaries belong in `results/`; full matrices and disposable runs do not.
The scripts contain no changes to either code's numerical implementation.
