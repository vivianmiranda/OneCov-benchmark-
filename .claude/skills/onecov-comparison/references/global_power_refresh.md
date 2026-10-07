# Global covariance power-table refresh

The adopted production preparation refines the original 1500-node logarithmic
wavenumber grid with a natural cubic spline, then keeps linear interpolation
in the C readers. At accuracy boost 1, `power_accuracyboost=8` gives
`(1500-1)*8+1 = 11993` nodes. This changes all installed linear, cb-linear
and nonlinear spectra, not only a four-halo diagnostic.

`initialize_forecast` (through each project's `initialize`) returns these
already installed tables in its existing dictionary format. Its resolved
settings record `power_refinement`. Do not call
`covariance.power.refine_power_tables` again on that return value.

## Existing-script audit

- `prepare_lsst_y1.py` and `complete_shear.py` already use the actual lengths
  of the returned k and z arrays. No script hardcodes a 1500-node power grid.
  Occurrences of 1500 in Gaussian case names refer to multipoles.
- The new export records all five power arrays, their shapes and float64
  byte hashes in `power_tables.npz` and its manifest. It requires 11993
  installed nodes for this refresh and records the production preparation
  source hash. The exact same exported nonlinear samples go to OneCov.
- `compare_halo.py`, `compare_trispectrum.py` and the complete CoCoA run
  record their installed power fingerprints too. The runner checks them
  against the shared export. Existing complete-matrix shared-Pmm checks
  remain intact; no grid-equality or tolerance check was removed.
- OneCov's supplied Pmm changes its nonlinear projection input. Its native
  halo/SSC **linear** spectrum still comes from its own hmf/CAMB preparation.
  This remains a native-model comparison, not identical linear power in
  both halo models. Shared-input formula/projection tests are labeled
  separately and use arrays captured from actual OneCov calculations.

## Sequential refresh

Run from this repository in the documented active Cocoa environment, after
the production API and targeted tests pass. The other interpreter is the
existing OneCov private environment. No installation or compilation occurs.

**Step 1:** inspect the commands. This performs no numerical calculation and
does not create the requested output directory.

```bash
python scripts/refresh_global_power.py --output work/global_power_11993_20261006
```

**Step 2:** execute the same sequence only after the numerical slot is free.

```bash
python scripts/refresh_global_power.py --output work/global_power_11993_20261006 --execute
```

The runner inherits `OMP_NUM_THREADS` from the shell (1–8 on this laptop),
requires `OPENBLAS_NUM_THREADS=1`, and uses one child at a time. It removes
Cocoa's Python import paths from the OneCov children's environment. Use
`--cocoa-python` or `--onecov-python` only to select an existing interpreter.

The 32 stages cover three Gaussian cases, SSC response and radial projection
(601 and 1201 nodes), halo ingredients, separated 1h–4h terms, shared cNG
projection, and complete small Fourier/real matrices. Every entry of G,
SSC, cNG and total is retained; the total positivity and generalized-mode
checks run with their existing plotting script. Both xi signs and their
cross-covariance appear in the real-space pilot.

Most stages have a ten-minute inner limit; complete Fourier and real runs
have fifteen- and thirty-minute limits. The runner has a two-hour overall
budget and starts a stage only when its entire deadline fits. It stops on
failure, preserves the log and folder, and records status in `refresh.json`.
Shell success is not a claim of general scientific convergence: inspect the
component differences, plot rendering and convergence diagnostics before
updating the README. Do not restart the whole sequence after a late failure;
use its saved command list to resume only unfinished stages in new folders.

The earlier OneCov real pilot took about 534 seconds including setup; the
Fourier pilot took about 82 seconds. These old values help budget the run,
but the denser input's reading/interpolation cost has not been measured.
Allow roughly 20–40 minutes initially, with the stated limits as safeguards.
All new outputs are accuracy evidence only. Collect fresh performance
measurements separately on a quiet machine, then update timing tables.

## Refinement and evidence boundaries

The ingredient checks keep global accuracy boost at 1 and change only
integration levels 0 and 2. This preserves their shared 11993-node grid.
Global boost 2 normally makes the refinement 16 and the grid 23985 nodes;
it therefore needs a matching fresh shared export. An explicit
`power_accuracyboost=4` at boost 2 is a separate fixed-power-grid control
and must be labeled as such, not silently substituted into an old scan.

The runner writes only a new work root. Old results, figures and current
README numbers remain untouched until the new records have been reviewed.
The old exhaustive aggregation script requires its full set of corner,
redshift, mass and cutoff studies; it is not suitable for this bounded
refresh. Broader claims from those studies require checking whether the
changed global spectra affect them, then rerunning the relevant controls.

## README follow-ups and timing

After the main sequence completes, the bounded follow-up refreshes the
sigma/abundance panels and the ten displayed connected-covariance controls:

```bash
python scripts/refresh_power_followups.py --baseline work/global_power_11993_20261006 --output work/global_power_11993_followups_20261006 --execute
```

It first checks that `diagnose_sigma.py` received exactly the installed
power fingerprint from the main input export. That script traces CAMB to
record the actual cosmology and amplitude, but its variance and sampled
power come from the initialized production C routines. The internal
sigma-table boost scan keeps the supplied 11993-node power fixed.

The connected controls preserve the old k, mass, redshift, radial and
corner settings one at a time. Their G+cNG mode ratios need refreshing:
even when native OneCov cNG is unchanged, G uses the new nonlinear Pmm.
Do not copy the old generalized-mode numbers into a new combined report.

Raw peak-height multiplicity/bias normalization identities are independent
of the power grid. Retain their evidence only after checking unchanged
source hashes, cosmology and fits; do not call old measurements new runs.

The native OneCov ingredient controls need more care. Its shear response
and higher-halo `Plin_spline` use native `mass_func.power`, not supplied
Pmm values, but `Setup.__consistency_checks_for_k_support_in_tabs` clamps
the hmf grid's support to the Pmm input. The shared minimum is
1.46184713e-5 h/Mpc, whereas the old standalone halo controls had no Pmm
and started at 1e-5. With 200 log samples, changing the lower endpoint
moves the entire grid. The cubic refinement preserves the original
shared-table endpoints: it did not itself create that range change.

Thus the current ingredients are made consistent with the complete
shared-Pmm pipeline. Cosmology, sigma8, survey parameters and n(z) hashes
were checked equal to the earlier complete export. The old native bias
divisor, corner, mass and SSC-control numbers must still be refreshed,
since those standalone configurations did not use the same native grid.
The follow-up now contains 40 stages, adding those controls. All its SSC
controls use 800 halo mass nodes; the earlier scan used 200. The full-fit
normalization identities, evaluated directly at supplied nu, are unchanged.

The prepared local command ledger is
`work/global_power_timing_plan_20261006.json`. It schedules 36 calls to
the existing timing scripts, without running them. It uses two fresh
processes per code for each Gaussian, SSC, cNG and halo-stage case, then
one complete Fourier and real-space pilot per code. First calls remain
separate from repeated means, and numerical checks run after measurement.

Use the new SSC601/1201 exports for the new timing rows, explicitly replacing
the old 300/2401-node rows rather than comparing unlike workloads. cNG keeps
300/601 nodes. Halo timing records the current installed-power fingerprint
and checks against the new integration-level-2 halo archive. Complete
timings omit `--accuracy-only`; no profiled export time enters these tables.
Run this ledger only after correctness and compiler jobs stop, sequentially,
with the documented eight OpenMP threads and single-thread BLAS setup.

Before timing, compare the committed preparation helper with the archived
inputs. The helper's file hash changed between the comparison export and
core commit `1d4428d`; do not infer numerical equivalence from that history.

```bash
python scripts/check_global_power_parity.py --native ../covariance_reference/four_halo_power_study/global_1/original_power.npz --dense work/global_power_11993_20261006/inputs --output work/global_power_11993_source_parity.json
```

This verifies the original and dense archives against their saved manifests,
runs the current helper on the original arrays, and requires exact float64
hashes and shapes for all five returned arrays. It records both source hashes
without rewriting either original manifest. No CAMB or covariance is called.
The timing results must retain their own installed-power fingerprints and
be checked against the same dense export, even after this parity check passes.

The stage timers record initialization and first-use calls separately from
repeated component calculations. The full CoCoA timer additionally records
its individual construction stages. OneCov's unmodified complete entry point
exposes setup and construction totals; those totals do not separately measure
halo-table generation and projection inside its calculation.

## Completed accuracy publication

Both runners completed: 32 main stages and 40 follow-up stages. Reviewed
records are in `results/global_power_20261006/`; `refresh.json` hashes every
published record and figure and identifies the unchanged supplied-nu CoCoA
normalization evidence separately from the newly run power-dependent tests.
Nine scientific figures were refreshed. No timing from these accuracy runs
is presented as a performance measurement; stale timing tables were removed
from the README pending the sequential timing run.

The archived-input parity check passed for all five arrays. Its record
preserves the original helper hash `43d4dc8c26d98f4b4727b26044fae182932bddf3635df5f75fa52b440f38243a`
and current hash `0456a460e4b5a48ea843836c8ccad4fb9b53e0a81771828335d6a3d65f0646b0`.
Exact installed-array equality, rather than assumed source-file equivalence,
connects the completed comparisons to committed core `1d4428d`.

The complete selected Fourier total has generalized variance ratios
0.9540178568–1.0069603903. The real-space symmetric totals give
0.9847882397–1.0234796978. Both are positive definite under their explicitly
stated eigenvalue conventions; these cross-code results are not Fisher or
full-survey convergence tests. The native OneCov diagonal angular-corner
and connected redshift refinements remain unresolved as documented in the
README. Project-wide production validation subsequently passed: 575 tests
across 139 modules, covariance-disabled checks and all seven notebooks.

## Completed quiet timing publication, 2026-10-07

The 36-stage sequential run completed without failed commands, after
all compiler, regression and notebook jobs had finished. The launch check
found 89.62% CPU idle, ordinary GUI activity and no numerical workers.
Eight OpenMP threads and one BLAS thread were used throughout. The actual
environment (`OMP_PROC_BIND=close`, unset `OMP_PLACES`/`OMP_DYNAMIC`) is
retained in the records; do not replace it with the platform recipe values
when describing these measurements.

`scripts/collect_power_timings.py` verified the installed power hashes,
all shared input files, repeated component outputs and the native settings.
Its first audit incorrectly required identical report key sets: accuracy
reports omit timings and instead contain an accuracy-only note. The
collector now permits those specific metadata differences, never changes
numerical checks, and retains both original reports. No covariance or
timing calculation was rerun.

All four complete timed NPZs are byte-identical to the accuracy archives.
Every G, SSC, cNG and total entry, coordinate and signal was also compared
with `np.array_equal`; all pass. OneCov's real-space antisymmetry is retained.
The two fresh processes per component give bitwise-identical outputs.
Shared-input cross-code residuals are below 4.21e-15 of the variance scale.

The published `results/global_power_20261006/timings.json` retains all
32 component records, four complete reports, original accuracy reports,
source hashes and command ledger. Component means pool 62 batch averages
(22 for halo moments), with sample standard deviation and separate process
means. Setup and first calls remain separate. Complete pilots are single
fresh-process measurements; never add warm component means to reconstruct
their construction times.

| Complete selected matrix | CoCoA construction | OneCov construction |
| --- | ---: | ---: |
| Fourier 100x100 | 20.9901 s | 71.5482 s |
| Real-space 16x16 | 45.7640 s | 491.4991 s |

Including setup gives 21.4833/90.6587 s for Fourier and
46.2743/516.7838 s for real space. Native halo moments take
37.59/2377.39 ms. In the supplied-trispectrum 8x8 projection OneCov is
faster: 0.0971 versus 0.2029 ms (300 nodes), 0.1065 versus 0.3336 ms
(601 nodes). Keep this result visible alongside CoCoA's faster SSC and
Gaussian component times. Native halo prescriptions and real-space
transforms still differ, and native cNG convergence remains unresolved.
These small cases do not measure a full 1560-entry OneCov forecast.

### Four-pilot cooling rerun, 2026-10-07

The user requested repeating only the four complete timings after setting
the fans to maximum and letting the laptop cool. That condition was
confirmed by the user, not measured by a thermal sensor. The prelaunch
check found 96.48% CPU idle and no other numerical/compiler jobs. The
same commands, scripts, eight-thread environment and inputs were used,
with fresh output folders; all earlier results remain intact.

| Complete pilot | CoCoA construction | OneCov construction |
| --- | ---: | ---: |
| Fourier 100x100 | 22.1845 s | 70.3927 s |
| Real-space 16x16 | 48.0055 s | 511.5998 s |

Changes from the earlier single runs are +5.69%/-1.61% for Fourier and
+4.90%/+4.09% for real space, in CoCoA/OneCov order. There is no systematic
speedup. Do not attribute these changes to temperature from one run per
condition, or select the fastest individual results from the two sets.
The README uses the complete new set. Repeated component timings were not
rerun and retain their own original records.

`scripts/collect_complete_timings.py` checks each matrix against both its
accuracy archive and earlier timing: every G, SSC, cNG, total, coordinate
and signal array is bitwise equal. Native settings and power fingerprints
also match. `results/global_power_20261006/complete_retiming_20261007.json`
preserves both sets of measurements and their verification. The missing
denominator brace in the README's matrix-difference equation was corrected
separately as f721f73; all ten display equations have balanced braces.

The subsequent full-survey discussion produced an unverified days-to-weeks
extrapolation from the native transform loops, 5-source tomography and
26 angular bins. It was withdrawn as a runtime estimate: proportional
cost per tomographic column has not been measured. Do not publish it as
a benchmark or revive it as an established prediction. CoCoA's shared
halo tables explain its inexpensive growth from a small pilot to full
LSST; OneCov's repeated real-space integrations need their own scaling
measurement. A bounded one/two-source and angular-bin study is the next
useful diagnostic, but was not part of the authorized four timing reruns.
