# TVSR 0.1.3

## Documentation hygiene

- Stripped BePe placeholder roxygen blocks (22 occurrences of
  `Convert a factor to numeric` / `fac2num(x)`) across `R/sdsr.r`,
  `R/lksr.r`, `R/spatialMatrix.r`, `R/plot.r`, `R/loglikelihood.r`,
  `R/complements.r`. Wrote minimum-correct `@title` / `@description` /
  `@param` / `@return` headers for the 11 exports that previously shared
  the placeholder (`SDSR` already had a proper docstring from 0.1.2);
  marked 8 internal helpers `@noRd` so they no longer generate stale
  `man/*.Rd`.
- Added `R/TVSR-package.R` declaring `@importFrom` for
  `stats::{df, median, nlm, optimize, rexp, rnorm, runif, sd}`,
  `utils::{modifyList, tail}`,
  `graphics::{grid, layout, legend, lines, par, plot.new}`, and
  `grDevices::{dev.off, pdf}`. Also calls `utils::globalVariables()` to
  silence "undefined global" NOTEs originating from dormant / orphan
  code paths in `R/2w.r`, `R/unused.r`, `R/test.r` pending their
  audit-driven removal (tracked in TVSR/TODO.md PR B).
- Added `R/data.R` with roxygen for the nine exported data objects:
  `df`, `dfPolicyRate`, `dfPolicyRateShadow`, `W`, `w`, `Y`, `time`,
  `sumTrade`, `Cftvw`. `w` shares an Rd topic with `W` via `@rdname` to
  avoid the `w.Rd` / `W.Rd` case-insensitive filesystem collision.
- Regenerated `man/` and `NAMESPACE` from scratch; the 13 stale
  internal / orphan Rd files (`h`, `hinv`, `K`, `Z`, `loglik`,
  `loglikStatic`, `loglikTVRhoCond`, `loglikTVRhoVarTrd`, `pause`,
  `SRllstatic`, `SRllTV`, `SRllTVCond`, `SRlocal`) were removed.
- Fixed non-ASCII characters in R source and data: French inline
  comments in `R/sdsr.r` translated to English; the degree-sign glyph
  `°` in three `cat()` strings removed; `"Côte d'Ivoire"` in
  `data/dbMonetaryPolicyInterestRates.rdata` re-encoded Latin-1 →
  UTF-8.
- Rewrote the DESCRIPTION `Description:` field to a formal
  CRAN-compliant wording and added the SSRN reference
  `<https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6589023>`.
  Synced the manual `Author:` / `Maintainer:` fields to match
  `Authors@R`-derived values.

## Bugs fixed

- `R/2w.r`: the dormant `loglikStaticAll2W` was calling
  `loglikStaticAll(..., omegaRho = ..., df = ...)` with two arguments
  that no longer exist on the current signature. Patched to a no-op
  stub matching the live signature; the function is not wired into any
  export and needs a proper port from p2 to become functional
  (tracked in TVSR/TODO.md P3).

## CI

- `.github/workflows/R-CMD-check.yaml`: tightened `error-on` from
  `"error"` to `"warning"`. `--no-examples` is still set; PR B (README
  + examples) will remove it once every export has a runnable
  `\examples{}` block.

## Status

- Local `R CMD check --no-examples --no-manual`: **0 ERROR, 0 WARNING,
  0 NOTE** on Linux release.

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
