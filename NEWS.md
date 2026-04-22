# TVSR 0.1.2

## Backwards-compatible additions

- `SDSR()` now exposes a four-spec decomposition through its `model`
  argument, which has been repurposed from the vestigial string
  `model = "trend"` into a named list:
  `model = list(trd, var)`, each field in `c("tv", "const", "gas")`.
  - `"tv"` (per-period analytical profile) is the 0.1.1 default; the
    `(tv, tv)` combination is bit-identical to 0.1.1.
  - `"const"` fits a pooled scalar nuisance by profile likelihood given
    the rho path: closed-form for `(tv, const)` and `(const, const)`;
    1-D Brent over the concentrated likelihood for `(const, tv)`.
  - `"gas"` (score-driven nuisance dynamics) is reserved for a future
    release and currently errors with an informative message.
- Two new internal helpers in `R/sdsr.r` underpin the four specs:
  `loglikTVRhoCond_mixed` (single-W) and `loglikTVRhoCond_mixed_tvW`
  (list-of-W). The originals `loglikTVRhoCond` and `loglikTVRhoCondtvW`
  are kept in place unchanged as reference implementations.
- First `tests/testthat/` suite and `.github/workflows/R-CMD-check.yaml`:
  tests cover `(tv, tv)` bit-identity, well-formed output on all four
  specs, flat paths for `"const"` specs, and the `"gas"` error. CI runs
  `R CMD check` on ubuntu (release + oldrel-1), macos, and windows.

## API changes

- `SDSR()`'s `model` argument changed type from character string
  (`"trend"`, never read in 0.1.1) to named list. Grep confirmed no
  in-tree or p2 caller passes `model = ...`, so this is a silent
  breaking change. Calling `SDSR(..., model = "trend")` now errors at
  `match.arg()` with a message pointing at the new list form.

## Bugs fixed

- `loglikStatic()` in `R/loglikelihood.r` previously hardcoded
  `KK = 1` (scalar) when `kernel = "uniform"`, which broke broadcasting
  for multi-period calls (`Nt > 1`). Now `KK = rep(1, Nt)`. All existing
  TVSR callers pass `Nt = 1` and are unaffected; the fix unblocks the
  p8 `srstatic_tvrho` estimator (plug-in nuisance + local TV rho over a
  uniform-kernel window) that would otherwise peg to the optimiser
  boundary.

## Documentation

- `SDSR()` gained a proper roxygen docstring describing the model, the
  four `model = list(trd, var)` spec options, and the profile-likelihood
  approach (replacing the generic stub that shipped in 0.1.1).
- A `FIXME` comment above the internal `SRllTV` documents the known
  `sample$omegaRho = ...` crash (root cause: `sample` is never
  initialised as a list) and the parallel bug in `p2`. `SRllTV` remains
  `@noRd` and dormant pending p2 archaeology in a future PR.

## Notes

- Vestigial `optim` argument on both `SDSR()` and `SRllTV` retained for
  signature compatibility; flagged for deprecation in a future release.
- Numerical verification: the numerical-equivalence harness at
  `p8_spatial_regression_bias/tvsr_numerical_equivalence_tests/` reports
  `max|diff| = 0` for the default SDSR, LKSR, and SRstatic paths.

# TVSR 0.1.1

## Backwards-compatible additions

- `SDSR()`, `LKSR()`, `SRstatic()`, and `crossValidLKSR()` now accept an
  `mc.cores` argument controlling the number of forked workers used in their
  internal grid searches. Default: `max(1L, parallel::detectCores() - 2L)`,
  which leaves 2 physical cores free for the rest of the workstation.
- Previously these functions hardcoded `parallel::detectCores()` internally.
  The new behaviour is a strict superset: callers that do not pass
  `mc.cores` still get a sensible default, and numerical outputs are
  unchanged (grid evaluations are deterministic, so worker count affects
  only wall-clock time). Verified by the numerical-equivalence test in
  `p8_spatial_regression_bias/tvsr_numerical_equivalence_tests/`.
- `LKSR()` and `crossValidLKSR()` now pass `mc.cores = 1L` to their internal
  `SRstatic()` calls to prevent nested parallelism when the outer loop is
  already forked.
