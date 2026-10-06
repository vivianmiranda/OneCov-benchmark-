---
name: onecov-comparison
description: Study and run OneCovariance against CoCoA for the LSST Y1 covariance benchmark, including matched inputs, bounded component tests, convergence and timing. Use for work in OneCov-benchmark-, not for changing production CosmoLike kernels.
---

# CoCoA–OneCovariance comparison

LSST Y1 is the benchmark. Begin with one source bin in Fourier space;
extend to the two-lens pilot and cross blocks after that calculation works.
The aim is to explain agreement and differences, not to make either code
win by choosing a cheaper or different physical model. **Do not optimize
OneCovariance.** Establish accuracy first, then time the implementations
at their validated settings. Source-level performance observations only
help budget affordable tests; they are not a request to rewrite the code.
Do not put recommendations for speeding up OneCovariance in the README.
Human documentation should explain matched inputs and relative differences.

Read [the source study](references/onecov.md) before changing input
formats, bias choices, numerical controls or component extraction. Use
the repository README for the current runnable commands and milestones.
Study the local OneCovariance checkout; do not fetch another implementation
or substitute online descriptions for the code actually being tested.

## Keep the comparison controlled

- Record revisions, input hashes, package versions and resolved settings.
  A configuration name alone is insufficient provenance.
- Begin with massless neutrinos, Limber, zero IA and matching survey area,
  source/lens distributions, number densities and per-component shape noise.
  Export the LSST project inputs rather than recreating its distributions.
- Use OneCovariance's supplied bias table for linear-bias galaxy tracers.
  Its default HOD model is a different calculation. Check that the table
  was actually consumed; the reader can otherwise fall back to HOD.
- First isolate Gaussian assembly with shared angular spectra. A native
  spectrum run is a separate test, since projection and matter-power
  differences enter there. Preserve both labels in results.
- Keep G, SSC and cNG separate. Set `split_gauss = True`; otherwise OneCov
  can fold SSC into cNG. Do not infer component agreement from the total.
- Match halo mass definitions, mass functions, profiles, concentration,
  damping and mass cutoffs before interpreting cNG differences as errors.
  The separate 2h/3h/4h export is not yet a public OneCov API.
- Never silently modify or monkey-patch either code to improve agreement.
  Propose a small, explicitly labeled diagnostic adapter when necessary.

## Laptop execution and accuracy

- Follow the main Cocoa and LSST Y1 READMEs for Cocoa activation,
  covariance compilation and platform-specific OpenMP/BLAS settings.
  Source `start_cocoa.sh` from `cocoa/Cocoa` in the Cocoa conda environment
  before returning to the benchmark directory. Do not replace this setup
  with a generic list of thread variables or add Apple vecLib settings
  as a Cocoa requirement.
- Installation follows Cocoa's schema: pinned choices in
  `set_installation_options.sh`, sourced setup/compile/start/stop scripts,
  the Cocoa Conda base and a repository-private `.local` environment.
  Downloads belong to setup; compilation is offline. Keep the two runtime
  environments in separate terminals. Do not invent a different setup
  scheme for each comparison repository.
- Use Python runners for numerical comparisons. Run one numerical job
  at a time. Take timings on a quiet machine, with BLAS fixed to one and
  at most eight CPU workers here. Derive worker count from
  `OMP_NUM_THREADS`; do not put an independent thread count in survey YAML.
- Start with a wall-time limit. A small output matrix may still initialize
  large shared tables: measure setup separately from assembly before
  predicting a full-survey runtime. Preserve failed and timed-out logs.
- Do not run a full OneCov matrix or a broad scan by default. Refine one
  selected block or ingredient, then decide whether expansion is affordable.
- Before broader agreement claims, cover low and high wavenumbers, more
  than one redshift, and off-diagonal entries. Keep this execution guidance
  here; the README should concentrate on results and their reproduction.
- Small real-space outputs still need converged multipole integration.
  Keep Cocoa's physical high-ell range when testing its real-space kernels;
  a cutoff suitable for a Fourier pilot does not certify real space.
- Check each code against its own finer settings. Compare diagonal
  variances, correlations and generalized variance eigenvalues, with the
  same ordering and cuts. Require positive total matrices; do not require
  each connected component to be positive definite by itself.
- Fisher convergence requires named parameters and derivatives. The
  data-vector delta-chi-squared criterion is not a covariance criterion.
- Full-construction timings include first-use shared tables. Reusing them
  is a separately labeled projection timing. No MCMC-style cache warmup
  may be hidden in the headline full-construction measurement.

## Changes and documentation

Keep scripts small and explicit. Document units, array axes, bin numbering
and what each major loop computes physically. Separate stages with blank
lines and explanatory comments. Avoid recovery frameworks for unusual
inputs when a clear guard suffices.

After each completed component, review the code for scientific correctness
and student readability before adding the next one. Test the actual parser
and smallest numerical case when dependencies are available; state clearly
when only preparation or parsing was tested.

Keep the human README self-contained; do not send readers to this skill
or a private study directory. Follow Cocoa's numbered Step format: give
each command its own step and code block. A continued command may span
lines; do not combine independent commands or semicolon chains in a box.
Commit coherent progress locally. Never push.
Installing missing packages requires explicit authorization; use an
isolated environment and do not upgrade the Cocoa environment.
