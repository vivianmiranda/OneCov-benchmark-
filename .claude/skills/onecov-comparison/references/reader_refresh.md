# Log-domain reader comparison refresh — 2026-10-08

## Scope and provenance

CoCoA adopted a log-domain linear-power reader for the connected
covariance (shared log10 k grid plus one per-shell shift), a per-call
z slice of the bilinear table read, and a blocked fusion that feeds the
tree averages directly from that reader. The merged core is `23b0127`;
the lsst_y1 project is `f7168af`; interface SHA-256
`3c361f80d018377cda141fd710e75adf27befe89ab945d371bd23c26cb27f1fb`.
Against the pre-adoption archive, Gaussian, SSC and signal are bitwise
identical and cNG entries move at most 4.3e-19 absolute, so only the
complete-matrix comparisons and every timing claim needed refreshing.
The ingredient comparisons and the component ms-scale timing tables
from 2026-10-06/07 are retained with their original record paths.

The refreshed artifacts live in `results/reader_refresh_20261008/`
(`fourier_shear.json`, `real_shear.json`, `pilot_timing_20261008.json`,
`full_lsst_timing_20261008.json`) with runs preserved under
`work/reader_refresh_20261008/`, including the full-LSST CLI logs in
`full_lsst_logs/`. The unchanged OneCov revision is
`311c2cfbe9d584d9d29abd3d8edaeea673470b7f`.

## The fresh export, and why two CoCoA run sets exist

The archived 2026-10-06 OneCov outputs failed their provenance check
against a fresh export: the exported Pmm grid shape changed (160500 to
155909 rows) after the project's dyadic z-grid adoption, with CAMB pins
identical. OneCov was therefore rerun from the fresh export
`work/lsst_y1_complete_20261008` in both spaces. Never reuse an
archived OneCov run whose `input_manifest_sha256` does not match the
current export; the comparison script enforces this.

Two CoCoA pilot run sets exist per space, and both are preserved:

- `{space}_cocoa` — timed runs against the archived shared inputs
  (`work/global_power_11993_20261006/inputs`), launched first on a quiet
  sequential machine. These provide the README timing table entries.
- `{space}_cocoa_v2` — comparison runs from the fresh export, published
  in the comparison records.

Every saved component of the timed runs is `array_equal` to its
comparison run (the CoCoA pilot reads the project's installed power, not
the exported files), so the quoted timings describe exactly the
published matrices. The verification is re-executed and embedded in
`pilot_timing_20261008.json` by the session's builder script.

## Results

- Fourier 100 by 100: total entry metric 0.328250% (was 0.328315%);
  generalized variance ratios 0.954018–1.006960; both totals positive
  definite. Real 16 by 16: total 2.125377% (was 2.124543%); ratios
  0.984788–1.023480; OneCov's largest antisymmetry 0.00114% (Gaussian).
  All README accuracy tables are unchanged at displayed precision; only
  record links and dates moved.
- Pilot timings (one fresh process, eight OpenMP threads, wall power,
  sequential): Fourier CoCoA 0.51/12.81 s versus OneCov 19.36/71.34 s
  (construction 5.57x, combined 6.81x). Real CoCoA 0.53/23.54 s versus
  OneCov 26.17/493.15 s (20.95x, 21.58x). The fan state was not
  instrumentally confirmed for these runs, so the README now says
  "measured sequentially" instead of repeating the fans-at-maximum
  wording of the 2026-10-07 retiming.
- Full LSST Y1 (1560 by 1560, project CLI example): 30.68/28.38/29.34 s
  construction, mean 29.47 s, initialization about 0.5 s. Small-to-full
  factor 29.47/23.54 = 1.25.

## Housekeeping

The never-referenced `wynn40_vs_native6*` and `wynn40_vs_production4*`
JSON/figure files sitting untracked in `results/` were moved to
`work/wynn40_orphaned_20261008/`: they were disposable intermediates of
the cutoff study, superseded by the committed `wynn_stability_vs_*`
records, and nothing in the README or the Wynn checkpoint cites them.
