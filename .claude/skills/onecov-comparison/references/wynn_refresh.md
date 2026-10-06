# Current-Wynn comparison refresh — 2026-10-06

## Scope and provenance

The current numerical implementation is core `d95867f`: FFTLog bias 0.8,
variance domain through k=1e25 h/Mpc, covariance mass integrals through
1e-40 Msun/h, and guarded Wynn acceleration of I11. The later core hashes
`26a0f97` and `2fb2821` in the saved reports include documentation changes.
Every refreshed CoCoA comparison uses interface SHA-256
`22e2e3499ac4de7bdcd52eaee62fb148a30d3773d1c49e684d827e8943d079cc`.

All 38 stages in `work/wynn_refresh_20261006/status.txt` exited zero. These
were accuracy-only runs at two OpenMP threads alongside the separate
all-project regression; their elapsed times are not benchmarks. The
all-project regression remains a separate acceptance gate.

The unchanged OneCov revision is
`311c2cfbe9d584d9d29abd3d8edaeea673470b7f`. Twenty-nine saved native metadata
hashes were checked before reusing those results. Shared-input tests reuse
the same spectra, responses and trispectra; native ingredient and complete
CoCoA shear matrices were recomputed. `results/wynn_refresh_20261006.json`
records native metadata hashes, refreshed report hashes and all nine
PNG/PDF figure pairs. The local runner and collector preserve their own
source fingerprints in the saved reports.

Previous results were preserved in
`work/wynn_refresh_20261006/previous_records/`, including
`before_aggregation/`. Earlier committed records remain in Git history.
Do not put historical cutoff studies back in the public README.

## Results and limits

- Shared Gaussian: largest variance-mode differences are 0.000017%,
  0.000465% and 0.000046%, consistent with native OneCov text precision.
- Shared SSC: maximum variance-scaled residual is 6.47e-15 across six
  configurations. Shared cNG projection remains at floating-point precision.
- Native Fourier shear (100 by 100): the total's largest entry difference,
  normalized by OneCov's total diagonal rms product, is 0.328315%. Total
  generalized variance ratios range from 0.954008 to 1.006961.
- Native real shear (16 by 16, xi+/xi- including their cross): the total
  entry metric is 2.124543%; generalized ratios span 0.984781–1.023477.
  The original matrices are preserved; only eigenvalue diagnostics use
  symmetric parts. Both totals' symmetric parts are positive definite.
- Native sigma differences are 0.03230%, 0.03036% and 0.02950%. The
  abundance normalization difference still grows toward z=1; Wynn does
  not change the fitted multiplicity amplitude or enforce both sum rules.
- CoCoA halo level 0 to 2 changes I11/I02/I12/response by at most 2.80e-7
  fractionally over 0.001–10 h/Mpc. Upper mass-panel nodes are 96 to 256;
  low-mass-panel nodes are 32 to 128. The separate native trispectrum check
  changes all terms by less than 3.92e-6 fractionally.
- OneCov's unequal-k 2h partition discrepancy and unresolved diagonal
  angular-cutoff convergence remain. No numerical source was modified by
  this refresh; it does not certify full-survey or Fisher convergence.

The sigma diagnostic's dense filter still integrates 1e-7–1e5 h/Mpc,
while production FFTLog extends to 1e25. The README now distinguishes
these domains; this diagnostic samples only resolved masses 1e10–1e15.
Do not present it as a tiny-mass-tail validation.

## Remaining timing work

The README retains explicitly dated supplied-table projection timings,
which do not generate native halo quantities. The old CoCoA native-halo
and full real-space timing claims were removed from the current tables;
those entries are pending a quiet sequential measurement with the current
mass panels and Wynn implementation. Existing timing JSON records are
preserved, and native OneCov timings remain identified by their date.
Do not substitute the concurrent refresh's elapsed times.

## Documentation checks

All nine PNG figures were inspected visually. Refreshed reports contain
finite JSON numbers; report/figure hashes and local links were checked.
The README parses as 25 Markdown tables, with 15 explicit anchors and 36
local file links resolving. No independent-reference claims or historical
cutoff studies were added to its comparison narrative.
